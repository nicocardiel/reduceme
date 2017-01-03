C------------------------------------------------------------------------------
C Version 26-April-2000                                         file: findarc.f
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
C Program: findarc
C Classification: wavelengths
C Description: Interactive arc line identification.
C
Comment
C
C Busca lineas en un espectro
C
        PROGRAM FINDARC
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUEBEG,TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXL
C NMAXL tambien en subrutinas FINDLINE, SUBPLOT1,2,3, MARCALINEA, FINDLINEA
        PARAMETER(NMAXL=5000)
C
        INTEGER NIDEN,NSAVE
        INTEGER I,K,I0,I1,I2,IC,IMIN,LRED
        INTEGER NL,NDEG,NDEG0
        INTEGER N0,N1,N2,N3
        INTEGER NWX,NWY
        INTEGER NB,NSIDE
        INTEGER NDIV,NC
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER L1,L2
        REAL XARC(NMAXL),LAMBDA(NMAXL)
        REAL XARC0(NMAXL),LAMBDA0(NMAXL)
        REAL LAMBDATABLE(NMAXL),LMIN
        REAL SX(NCMAX),SY(NCMAX)
        REAL XC,YC
        REAL STWV0,DISP0,POL
        REAL SIGMAY(NMAXL),CHISQR,A(20),B(20)
        REAL YMIN,YMAX
        REAL X0,FINDMAX
        CHARACTER*1 COPC,CH,CN0
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,OUTFILE,ARCFILE,REDUCEMEDIR
        CHARACTER*80 CLAMBDATABLE(NMAXL)
        LOGICAL FITDONE
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LBEXIST
C
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKS/SX,SY
        COMMON/BLKNSIDE/NSIDE
        COMMON/BLKFIND1/NL
        COMMON/BLKFIND2/LAMBDATABLE
        COMMON/BLKYMMM/YMIN,YMAX
        COMMON/BLKPL1/A,B
        COMMON/BLKPL2/NIDEN,NDEG
        COMMON/BLKPL3/XARC,LAMBDA
        COMMON/BLKINFILE/INFILE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        THISPROGRAM='findarc'
        CALL WELCOME('6-December-1996')
C
        NDEG=1
C
        FITDONE=.FALSE.
        DO I=1,NMAXL
          SIGMAY(I)=0.
        END DO
C
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
C
        WRITE(*,100)'Spectrum file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(TRUELEN(OBJECT).GT.0)THEN
          INFILE=INFILE(1:TRUELEN(INFILE))//' ['//
     +     OBJECT(1:TRUELEN(OBJECT))//']'
        END IF
        READ(20) (SY(I),I=1,NCHAN)
        CLOSE(20)
C
        DO I=1,NCHAN
          SX(I)=REAL(I)
        END DO
C
        CALL GETENV('reduceme_dir',REDUCEMEDIR)
        LRED=TRUELEN(REDUCEMEDIR)
C
        WRITE(*,100)'Enter file with wavelenghts for current arc type '
        ARCFILE(1:LRED)=REDUCEMEDIR(1:LRED)
        ARCFILE(LRED+1:)='/files/cuar.ldat'
        ARCFILE=INFILEX(25,ARCFILE,0,0,.0,.0,3,.FALSE.)
        I=0
8       CONTINUE
        IF(I.EQ.NMAXL)THEN
          WRITE(*,101)'ERROR: No. of lines too large.'
          CLOSE(25)
          STOP
        END IF
        READ(25,101,END=9) CLAMBDATABLE(I+1)
        IF(CLAMBDATABLE(I+1)(1:1).NE.'#')THEN
          READ(CLAMBDATABLE(I+1),*,END=9) LAMBDATABLE(I+1)
          I=I+1
        END IF
        GOTO 8
9       CONTINUE
        NL=I
        CLOSE(25)
        WRITE(*,110)'No. of line positions read: ',NL
C
        WRITE(*,100)'No. of channels at each side to fit maximum '
        NSIDE=READILIM('2',1,9999)
C
        CALL BUTTON(1,'[Z]oom',0)
        CALL BUTTON(2,'[W]hole',0)
        CALL BUTTON(2,'[W]hole',3)
        CALL BUTTON(3,'[P]ol.deg.',0)
        CALL BUTTON(4,'[S]ave',0)
        CALL BUTTON(4,'[S]ave',3)
        CALL BUTTON(5,'[L]oad',0)
        CALL BUTTON(6,'[Q]UIT',0)
        I1=1
        I2=NCHAN
        CALL SUBPLOT1(I1,I2)
C
        NIDEN=0
C
10      CONTINUE
        CALL RPGBAND(7,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        IF((CH.EQ.'Z').OR.(CH.EQ.'z'))THEN
          NB=1
        ELSEIF((CH.EQ.'W').OR.(CH.EQ.'w'))THEN
          CALL BUTTQEX(2,LBEXIST)
          IF(LBEXIST) NB=2
        ELSEIF((CH.EQ.'P').OR.(CH.EQ.'p'))THEN
          NB=3
        ELSEIF((CH.EQ.'S').OR.(CH.EQ.'s'))THEN
          CALL BUTTQEX(4,LBEXIST)
          IF(LBEXIST) NB=4
        ELSEIF((CH.EQ.'L').OR.(CH.EQ.'l'))THEN
          NB=5
        ELSEIF((CH.EQ.'Q').OR.(CH.EQ.'q'))THEN
          NB=6
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.0)THEN
          IF((YC.GT.YMAX).OR.(YC.LT.YMIN)) GOTO 10
          WRITE(*,110)'Next line to be identified is #',NIDEN+1
          I0=NINT(XC)
          IF((I0.LT.NSIDE+1).OR.(NCHAN-I0.LT.NSIDE+1))THEN
            WRITE(*,101)'WARNING: current channel is too near to '//
     +       'the spectrum edge.'
            GOTO 10
          END IF
          IF(NIDEN.GT.0)THEN
            DO I=1,NIDEN
              IF(INT(XARC(I)).EQ.I0)THEN
                WRITE(*,100)'WARNING: this line has been already '//
     +           'identified with wavelength: '
                WRITE(*,*)LAMBDA(I)
                WRITE(*,100)'Do you want to change this line '//
     +           '(y/n/d=delete-line) '
                COPC(1:1)=READC('n','ynd')
                IF(COPC.EQ.'y')THEN
                  WRITE(*,100)'New wavelength '
                  WRITE(CDUMMY,*)LAMBDA(I)
                  LAMBDA(I)=READF(CDUMMY)
                  GOTO 11
                ELSEIF(COPC.EQ.'n')THEN
                  GOTO 10
                ELSE
                  GOTO 34
                END IF
              END IF
            END DO
          END IF
          GOTO 35
34        IF(I.EQ.NIDEN)THEN
            NIDEN=NIDEN-1
          ELSE
            DO K=I,NIDEN-1
              XARC(K)=XARC(K+1)
              LAMBDA(K)=LAMBDA(K+1)
            END DO
            NIDEN=NIDEN-1
          END IF
          GOTO 11
35        DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PGMOVE(REAL(I0),YMAX)
            CALL PGDRAW(REAL(I0),YMAX-(YMAX-YMIN)/10.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,110)'Mouse at channel #',I0
          WRITE(*,100)'Is it right (y/n) '
          CH(1:1)=READC('y','yn')
          IF(CH.EQ.'n')THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(0)
              CALL PGMOVE(REAL(I0),YMAX)
              CALL PGDRAW(REAL(I0),YMAX-(YMAX-YMIN)/10.)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            GOTO 10
          END IF
          X0=FINDMAX(I0)
          WRITE(*,100)'Maximum at: '
          WRITE(*,*) X0
          IF(FITDONE)THEN
            WRITE(*,100)'Estimated wavelength: '
            POL=B(NDEG+1)
            DO I=NDEG,1,-1
              POL=POL*X0+B(I)
            END DO
            WRITE(*,*)POL
            CALL FINDLINE(POL,N1,N2,N3,N0)
            WRITE(*,101)'Nearest lines in file:'
            WRITE(*,100) '(1) lambda info, increment: '
            L1=TRUEBEG(CLAMBDATABLE(N1))
            L2=TRUELEN(CLAMBDATABLE(N1))
            WRITE(*,100) CLAMBDATABLE(N1)(L1:L2)
            WRITE(*,*) LAMBDATABLE(N1)-POL
            WRITE(*,100) '(2) lambda info, increment: '
            L1=TRUEBEG(CLAMBDATABLE(N2))
            L2=TRUELEN(CLAMBDATABLE(N2))
            WRITE(*,100) CLAMBDATABLE(N2)(L1:L2)
            WRITE(*,*) LAMBDATABLE(N2)-POL
            WRITE(*,100) '(3) lambda info, increment: '
            L1=TRUEBEG(CLAMBDATABLE(N3))
            L2=TRUELEN(CLAMBDATABLE(N3))
            WRITE(*,100) CLAMBDATABLE(N3)(L1:L2)
            WRITE(*,*) LAMBDATABLE(N3)-POL
            WRITE(*,100) 'Option (1/2/3/4=1+2/5=2+3/6=1+3/0=none/'
            WRITE(*,100) '9=insert-by-keyboard) '
            WRITE(CDUMMY,*)N0
            CN0(1:1)=READC(CDUMMY,'01234569')
            READ(CN0,*)N0
          ELSE
            N0=9
          END IF
          IF(N0.EQ.0)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(0)
              CALL PGMOVE(X0,YMAX)
              CALL PGDRAW(X0,YMAX-(YMAX-YMIN)/10.)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            GOTO 10
          END IF
          NIDEN=NIDEN+1
          CALL BUTTON(4,'[S]ave',0)
          IF(N0.EQ.9)THEN
            WRITE(*,100)'Wavelength at this channel'
            LAMBDA(NIDEN)=READF('@')
          ELSE IF(N0.EQ.1)THEN
            LAMBDA(NIDEN)=LAMBDATABLE(N1)
          ELSE IF(N0.EQ.2)THEN
            LAMBDA(NIDEN)=LAMBDATABLE(N2)
          ELSE IF(N0.EQ.3)THEN
            LAMBDA(NIDEN)=LAMBDATABLE(N3)
          ELSE IF(N0.EQ.4)THEN
            LAMBDA(NIDEN)=(LAMBDATABLE(N1)+LAMBDATABLE(N2))/2.
          ELSE IF(N0.EQ.5)THEN
            LAMBDA(NIDEN)=(LAMBDATABLE(N2)+LAMBDATABLE(N3))/2.
          ELSE IF(N0.EQ.6)THEN
            LAMBDA(NIDEN)=(LAMBDATABLE(N1)+LAMBDATABLE(N3))/2.
          END IF
          XARC(NIDEN)=X0
          WRITE(*,*)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL MARCALINEA(NIDEN,.FALSE.,.TRUE.,LCOLOR(ITERM))
          END DO
11        IF(NIDEN.GE.2)THEN
            WRITE(*,110)'> Polynomial degree...........: ',NDEG
            WRITE(*,110)'> No. of lines in fit.........: ',NIDEN
            CALL POLFIT(XARC,LAMBDA,SIGMAY,NIDEN,2,0,A,CHISQR)
            DISP0=A(2)
            STWV0=A(1)+DISP0
            WRITE(*,100)'> Mean dispersion (A/pixel)...: '
            WRITE(*,*)DISP0
            WRITE(*,100)'> Wavelength at channel #1....: '
            WRITE(*,*)STWV0
            WRITE(*,*)
            IF(NDEG.GT.1)THEN
              CALL POLFIT(XARC,LAMBDA,SIGMAY,NIDEN,NDEG+1,0,B,CHISQR)
            ELSE 
              B(1)=A(1)
              B(2)=A(2)
            END IF
            CALL SUBPLOT2
            FITDONE=.TRUE.
          END IF
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[Z]oom',5)
          WRITE(*,100)'Point #1: press mouse button...'
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          WRITE(*,101)' OK!'
          I1=NINT(XC)
          WRITE(*,100)'Point #2: press mouse button...'
          CALL RPGBAND(4,0,REAL(I1),0.,XC,YC,CH)
          WRITE(*,101)' OK!'
          I2=NINT(XC)
          IF(I2.LT.I1)THEN
            IC=I1
            I1=I2
            I2=IC
          END IF
          IF(I1.LT.1) I1=1
          IF(I2.GT.NCHAN) I2=NCHAN
          CALL SUBPLOT1(I1,I2)
          CALL BUTTON(1,'[Z]oom',0)
          CALL BUTTON(2,'[W]hole',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'[W]hole',5)
          I1=1
          I2=NCHAN
          CALL SUBPLOT1(I1,I2)
          CALL BUTTON(2,'[W]hole',0)
          CALL BUTTON(2,'[W]hole',3)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[P]ol.deg.',5)
          NDEG0=NDEG
          WRITE(*,100)'Polynomial degree '
          WRITE(CDUMMY,*)NDEG
          NDEG=READI(CDUMMY)
          IF(NDEG.GT.NIDEN-1)THEN
            WRITE(*,101)'ERROR: polynomial degree exceeds data number.'
            NDEG=NDEG0
            WRITE(*,110)'Polynomial degree forced to be: ',NDEG
          END IF
          IF(NDEG.GT.19)THEN
            WRITE(*,101)'ERROR: polynomial degree must be < 20.'
            NDEG=NDEG0
            WRITE(*,110)'Polynomial degree forced to be: ',NDEG
          END IF
          IF(NDEG.LT.1)THEN
            WRITE(*,101)'ERROR: polynomial degree must be > 1.'
            NDEG=NDEG0
            WRITE(*,110)'Polynomial degree forced to be: ',NDEG
          END IF
          WRITE(*,110)'> Polynomial degree...........: ',NDEG
          WRITE(*,110)'> No. of lines in fit.........: ',NIDEN
          CALL POLFIT(XARC,LAMBDA,SIGMAY,NIDEN,2,0,A,CHISQR)
          DISP0=A(2)
          STWV0=A(1)+DISP0
          WRITE(*,100)'> Mean dispersion (A/pixel)...: '
          WRITE(*,*)DISP0
          WRITE(*,100)'> Wavelength at channel #1....: '
          WRITE(*,*)STWV0
          WRITE(*,*)
          IF(NDEG.GT.1)THEN
            CALL POLFIT(XARC,LAMBDA,SIGMAY,NIDEN,NDEG+1,0,B,CHISQR)
          ELSE 
            B(1)=A(1)
            B(2)=A(2)
          END IF
          CALL SUBPLOT2
          CALL BUTTON(3,'[P]ol.deg.',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.4)THEN
C salvamos las lineas identificadas en orden de l.d.o.
          CALL BUTTON(4,'[S]ave',5)
ccc13        WRITE(*,100)'Output file name (save positions and '//
          WRITE(*,100)'Output file name (save positions and '//
     +     'wavelenghts) '
          OUTFILE=OUTFILEX(30,'fitlin.dat',0,0,0.,0.,3,.FALSE.)
          NSAVE=NIDEN
          DO I=1,NSAVE
            XARC0(I)=XARC(I)
            LAMBDA0(I)=LAMBDA(I)
          END DO
14        IMIN=1
          IF(NSAVE.GT.1)THEN
            LMIN=XARC0(1)
            DO I=2,NSAVE
              IF(XARC0(I).LT.LMIN)THEN
                LMIN=XARC0(I)
                IMIN=I
              END IF
            END DO
          END IF
          WRITE(30,*)NINT(XARC0(IMIN)),LAMBDA0(IMIN)
          IF(IMIN.LT.NSAVE)THEN
            DO I=IMIN,NSAVE-1
              XARC0(I)=XARC0(I+1)
              LAMBDA0(I)=LAMBDA0(I+1)
            END DO
          END IF
          NSAVE=NSAVE-1
          IF(NSAVE.GT.0) GOTO 14
          CLOSE(30)
          WRITE(*,110)'No. of lines saved: ',NIDEN
          CALL BUTTON(4,'[S]ave',0)
          CALL BUTTON(4,'[S]ave',3)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.5)THEN
          CALL BUTTON(5,'[L]oad',5)
ccc16        WRITE(*,100)'File with positions and wavelenghts '
          WRITE(*,100)'File with positions and wavelenghts '
          OUTFILE=INFILEX(32,'fitlin.dat',0,0,.0,.0,3,.FALSE.)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
          END DO
          I=0
17        CONTINUE
          READ(32,*,END=18)XARC(I+1),LAMBDA(I+1)
          I=I+1
          IF((XARC(I).GE.I1).AND.(XARC(I).LE.I2))THEN 
            IF((I1.EQ.1).AND.(I2.EQ.NCHAN))THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL MARCALINEA(I,.FALSE.,.FALSE.,LCOLOR(ITERM))
              END DO
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL MARCALINEA(I,.FALSE.,.TRUE.,LCOLOR(ITERM))
              END DO
            END IF
          END IF
          GOTO 17
18        CONTINUE
          CLOSE(32)
          NIDEN=I
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C
          IF(NIDEN.GE.2)THEN
            WRITE(*,110)'> Polynomial degree...........: ',NDEG
            WRITE(*,110)'> No. of lines in fit.........: ',NIDEN
            CALL POLFIT(XARC,LAMBDA,SIGMAY,NIDEN,2,0,A,CHISQR)
            DISP0=A(2)
            STWV0=A(1)+DISP0
            WRITE(*,100)'> Mean dispersion (A/pixel)...: '
            WRITE(*,*)DISP0
            WRITE(*,100)'> Wavelength at channel #1....: '
            WRITE(*,*)STWV0
            WRITE(*,*)
            IF(NDEG.GT.1)THEN
              CALL POLFIT(XARC,LAMBDA,SIGMAY,NIDEN,NDEG+1,0,B,CHISQR)
            ELSE 
              B(1)=A(1)
              B(2)=A(2)
            END IF
            CALL SUBPLOT2
            FITDONE=.TRUE.
          END IF
          CALL BUTTON(5,'[L]oad',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[Q]UIT',5)
          WRITE(*,100)'Do you really want to QUIT (y/n) '
          COPC(1:1)=READC('y','yn')
          IF(COPC.EQ.'y') GOTO 80
          CALL BUTTON(6,'[Q]UIT',0)
C------------------------------------------------------------------------------
        ELSE
          WRITE(*,101)'FATAL ERROR. Imposible option.'
          STOP
        END IF
C
        GOTO 10
C
80      CONTINUE
C
        CALL PGEND
C
        IF(NIDEN.EQ.0) GOTO 90
        WRITE(*,100) 'Do you want to see plots '//
     +   'with the identified lines (y/n) '
        COPC(1:1)=READC('n','yn')
        IF(COPC.EQ.'y')THEN
          WRITE(*,110) 'No. of channels: ',NCHAN
          WRITE(*,100) 'How many plots do you want to divide the '//
     +     'spectrum into '
          NDIV=READI('5')
          WRITE(*,100) 'No. of plots in X,Y '
          CALL READ2I('1,1',NWX,NWY)
C usamos un poco de solape en los plots
          CALL PGBEGIN(0,'?',NWX,NWY)
          NC=NINT(REAL(NCHAN)/REAL(NDIV))+1
          NC=NINT(REAL(NC)*1.2)
          DO I=1,NDIV
            IC=NINT((REAL(I)-0.5)/REAL(NDIV)*REAL(NCHAN))
            I1=IC-NC/2
            I2=IC+NC/2
            IF(I1.LT.1) I1=1
            IF(I2.GT.NCHAN) I2=NCHAN
            CALL SUBPLOT3(I1,I2)
          END DO
          CALL PGEND
        END IF
C
90      CONTINUE
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE SUBPLOT1(I1,I2)
        IMPLICIT NONE
        INTEGER I1,I2
        INCLUDE 'redlib.inc'
C
        INTEGER NMAXL
        PARAMETER(NMAXL=5000)
C
        INTEGER I
        INTEGER NIDEN,NDEG
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL SX(NCMAX),SY(NCMAX)
        REAL XARC(NMAXL),LAMBDA(NMAXL)
        CHARACTER*75 INFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKPL2/NIDEN,NDEG
        COMMON/BLKPL3/XARC,LAMBDA
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKS/SX,SY
        COMMON/BLKYMMM/YMIN,YMAX
        COMMON/BLKINFILE/INFILE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(I1.EQ.I2)THEN
          WRITE(*,101)'ERROR: invalid limits in subroutine SUBPLOT1.'
          WRITE(*,100)'(press RETURN to continue)'
          READ(*,*)
        END IF
        YMIN=SY(I1)
        YMAX=YMIN
        DO I=I1+1,I2
          IF(SY(I).LT.YMIN) YMIN=SY(I)
          IF(SY(I).GT.YMAX) YMAX=SY(I)
        END DO
        XMIN=REAL(I1)
        XMAX=REAL(I2)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(0.,1.,0.,0.55)
          CALL PGWINDOW(0.,1.,0.,1.)
          CALL PGSFS(1)
          CALL PGSCI(0)
          CALL PGRECT(0.,1.,0.,1.)
          CALL PGSCI(1)
          CALL PGVPORT(.1,.95,.1,.53)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
          CALL PGMTEXT('T',0.5,0.5,0.5,'file: '//INFILE)
          CALL PGIDEN_RED
          IF(LCOLOR(ITERM)) CALL PGSCI(3)
          CALL PGBIN(NCHAN,SX,SY,.TRUE.)
          DO I=1,NIDEN
            IF((INT(XARC(I)).GE.I1).AND.(INT(XARC(I)).LE.I2))THEN
              IF((I1.EQ.1).AND.(I2.EQ.NCHAN))THEN
                CALL MARCALINEA(I,.FALSE.,.FALSE.,LCOLOR(ITERM))
              ELSE
                CALL MARCALINEA(I,.FALSE.,.TRUE.,LCOLOR(ITERM))
              END IF
            END IF
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('channel','No. of counts',CHAR(32))
        END DO
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE SUBPLOT2
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER NMAXL
        PARAMETER(NMAXL=5000)
C
        INTEGER I,J
        INTEGER NIDEN,NDEG
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(20),B(20)
        REAL XARC(NMAXL),LAMBDA(NMAXL)
        REAL X(NMAXL),Y(NMAXL)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL POL
        REAL X1,X2,Y1,Y2
        REAL XV1,XV2,YV1,YV2
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKPL1/A,B
        COMMON/BLKPL2/NIDEN,NDEG
        COMMON/BLKPL3/XARC,LAMBDA
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        CALL PGQVP(0,XV1,XV2,YV1,YV2)
        CALL PGQWIN(X1,X2,Y1,Y2)
C
        YMIN=1.E6
        YMAX=-1.E6
        DO I=1,NIDEN
          POL=A(NDEG+1)
          DO J=NDEG,1,-1
            POL=POL*XARC(I)+A(J)
          END DO
          Y(I)=LAMBDA(I)-POL
          IF(Y(I).GT.YMAX) YMAX=Y(I)
          IF(Y(I).LT.YMIN) YMIN=Y(I)
        END DO
        IF(ABS(YMIN).LT.A(2)) YMIN=-A(2)
        IF(ABS(YMAX).LT.A(2)) YMAX=A(2)
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
        XMIN=1.
        XMAX=REAL(NCHAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGVPORT(.0,1.,0.57,0.85)
          CALL PGWINDOW(0.,1.,0.,1.)
          CALL PGSFS(1)
          CALL PGSCI(0)
          CALL PGRECT(0.,1.,0.,1.)
          CALL PGSCI(1)
C
          CALL PGVPORT(.1,.95,0.60,0.80)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('BCMST',0.0,0,'BCNST',0.0,0)
          CALL PGIDEN_RED
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGPOINT(NIDEN,XARC,Y,16)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        DO I=1,NMAXL
          X(I)=REAL(I-1)/REAL(NMAXL-1)*REAL(NCHAN-1)+1.
          IF(NDEG.EQ.1)THEN
            Y(I)=0.
          ELSE
            POL=B(NDEG+1)
            DO J=NDEG,1,-1
              POL=POL*X(I)+B(J)
            END DO
            Y(I)=POL-A(1)-A(2)*X(I)
          END IF
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(6)
          CALL PGLINE(NMAXL,X,Y)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGVPORT(XV1,XV2,YV1,YV2)
          CALL PGWINDOW(X1,X2,Y1,Y2)
        END DO
        END
C
C******************************************************************************
C
        SUBROUTINE SUBPLOT3(I1,I2)
        IMPLICIT NONE
        INTEGER I1,I2
        INCLUDE 'redlib.inc'
C
        INTEGER NMAXL
        PARAMETER(NMAXL=5000)
C
        INTEGER I
        INTEGER NIDEN,NDEG
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL SX(NCMAX),SY(NCMAX)
        REAL XARC(NMAXL),LAMBDA(NMAXL)
        CHARACTER*75 INFILE
C
        COMMON/BLKPL2/NIDEN,NDEG
        COMMON/BLKPL3/XARC,LAMBDA
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKS/SX,SY
        COMMON/BLKYMMM/YMIN,YMAX
        COMMON/BLKINFILE/INFILE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(I1.EQ.I2)THEN
          WRITE(*,101)'ERROR: invalid limits in subroutine SUBPLOT3.'
          WRITE(*,100)'(press RETURN to continue)'
          READ(*,*)
        END IF
        YMIN=SY(I1)
        YMAX=YMIN
        DO I=I1+1,I2
          IF(SY(I).LT.YMIN) YMIN=SY(I)
          IF(SY(I).GT.YMAX) YMAX=SY(I)
        END DO
        XMIN=REAL(I1)
        XMAX=REAL(I2)
        DX=XMAX-XMIN
        DY=YMAX-YMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
        CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
        CALL PGIDEN_RED
        CALL PGBIN(NCHAN,SX,SY,.TRUE.)
        DO I=1,NIDEN
          IF((XARC(I).GE.I1).AND.(XARC(I).LE.I2))THEN
            CALL MARCALINEA(I,.TRUE.,.TRUE.,.FALSE.)
          END IF
        END DO
        CALL PGLABEL('channel','No. of counts','file: '//INFILE)
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C*****************************************************************************
C
C LMODE=.TRUE. usamos todos los decimales disponibles en l.d.o.
C LMODE=.FALSE. redondeamos a las unidades en l.d.o.
        SUBROUTINE MARCALINEA(N,LMODE,LLDO,LLCOLOR)
        IMPLICIT NONE
        INTEGER N
        LOGICAL LMODE,LLDO,LLCOLOR
C
        INCLUDE 'redlib.inc'
        INTEGER TRUEBEG,TRUELEN
C
        INTEGER NMAXL
        PARAMETER(NMAXL=5000)
C
        INTEGER L1,L2
        REAL SX(NCMAX),SY(NCMAX)
        REAL YMIN,YMAX,DY
        REAL XARC(NMAXL),LAMBDA(NMAXL)
        REAL XX,YY
        CHARACTER*50 CDUMMY
C
        COMMON/BLKS/SX,SY
        COMMON/BLKYMMM/YMIN,YMAX
        COMMON/BLKPL3/XARC,LAMBDA
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DY=YMAX-YMIN
        XX=XARC(N)
        IF(LLCOLOR) CALL PGSCI(2)
        CALL PGMOVE(XX,YMAX)
        CALL PGDRAW(XX,YMAX-DY/10.)
        IF(LLCOLOR) CALL PGSCI(1)
        IF(LLDO)THEN
          YY=SY(INT(XX))+DY/25.
          IF(LMODE)THEN
            WRITE(CDUMMY,'(F8.1)')LAMBDA(N)
          ELSE
            WRITE(CDUMMY,'(I5)')NINT(LAMBDA(N))
          END IF
          L1=TRUEBEG(CDUMMY)
          L2=TRUELEN(CDUMMY)
          CALL PGSCH(.7)
          CALL PGPTEXT(XX,YY,0.,.5,CDUMMY(L1:L2))
          CALL PGSCH(1.)
          CALL PGMOVE(XX,YY-DY/100.)
          CALL PGDRAW(XX,SY(INT(XX))+DY/100.)
        END IF
        IF(LLCOLOR) CALL PGSCI(1)
        END
C
C*****************************************************************************
C
C Determina las tres lineas mas proximas (en el fichero con lineas del arco
C identificadas) a la estimacion X, realizada por un ajuste lineal.
C N1, N2, N3 son el numero de las 3 lineas mas proximas a X en la matriz
C LAMBDATABLE(). N0 es el numero de orden de la linea mas proxima (que puede
C ser cualquiera entre N1, N2 y N3).
        SUBROUTINE FINDLINE(X,N1,N2,N3,N0)
        IMPLICIT NONE
        INTEGER N1,N2,N3,N0
        REAL X
C
        INTEGER NMAXL
        PARAMETER(NMAXL=5000)
C
        INTEGER I
        INTEGER NL
        INTEGER NN1,NN2,NN3
        REAL LAMBDATABLE(NMAXL)
        REAL DX,DXMIN
        LOGICAL LNN1,LNN2,LNN3
C
        COMMON/BLKFIND1/NL
        COMMON/BLKFIND2/LAMBDATABLE
C------------------------------------------------------------------------------
C
C buscamos la linea mas proxima
        NN1=0
        DXMIN=1.E16
        DO I=1,NL
          DX=ABS(LAMBDATABLE(I)-X)
          IF(DX.LT.DXMIN)THEN
            DXMIN=DX
            NN1=I
          END IF
        END DO
        IF(NN1.EQ.0)THEN
          WRITE(*,101)'FATAL ERROR #1 in subroutine FINLINE.'
          STOP
        END IF
C buscamos otra vez la linea mas proxima (anulando NN1)
        NN2=0
        DXMIN=1.E16
        DO I=1,NL
          IF(I.NE.NN1)THEN
            DX=ABS(LAMBDATABLE(I)-X)
            IF(DX.LT.DXMIN)THEN
              DXMIN=DX
              NN2=I
            END IF
          END IF
        END DO
        IF(NN2.EQ.0)THEN
          WRITE(*,101)'FATAL ERROR #2 in subroutine FINLINE.'
          STOP
        END IF
C buscamos otra vez la linea mas proxima (anulando NN1 y NN2)
        NN3=0
        DXMIN=1.E16
        DO I=1,NL
          IF((I.NE.NN1).AND.(I.NE.NN2))THEN
            DX=ABS(LAMBDATABLE(I)-X)
            IF(DX.LT.DXMIN)THEN
              DXMIN=DX
              NN3=I
            END IF
          END IF
        END DO
        IF(NN3.EQ.0)THEN
          WRITE(*,101)'FATAL ERROR #3 in subroutine FINLINE.'
          STOP
        END IF
C ordenamos NN1, NN2 y NN3 (lo hacemos a lo bruto)
        LNN1=.TRUE.
        LNN2=.TRUE.
        LNN3=.TRUE.
C
        DO I=1,NL
          IF((I.EQ.NN1).AND.LNN1)THEN
            N1=NN1
            LNN1=.FALSE.
            N0=1
            GOTO 10
          END IF
          IF((I.EQ.NN2).AND.LNN2)THEN
            N1=NN2
            LNN2=.FALSE.
            GOTO 10
          END IF
          IF((I.EQ.NN3).AND.LNN3)THEN
            N1=NN3
            LNN3=.FALSE.
            GOTO 10
          END IF
        END DO
C
10      DO I=1,NL
          IF((I.EQ.NN1).AND.LNN1)THEN
            N2=NN1
            LNN1=.FALSE.
            N0=2
            GOTO 20
          END IF
          IF((I.EQ.NN2).AND.LNN2)THEN
            N2=NN2
            LNN2=.FALSE.
            GOTO 20
          END IF
          IF((I.EQ.NN3).AND.LNN3)THEN
            N2=NN3
            LNN3=.FALSE.
            GOTO 20
          END IF
        END DO
C
20      DO I=1,NL
          IF((I.EQ.NN1).AND.LNN1)THEN
            N3=NN1
            LNN1=.FALSE.
            N0=3
            GOTO 30
          END IF
          IF((I.EQ.NN2).AND.LNN2)THEN
            N3=NN2
            LNN2=.FALSE.
            GOTO 30
          END IF
          IF((I.EQ.NN3).AND.LNN3)THEN
            N3=NN3
            LNN3=.FALSE.
            GOTO 30
          END IF
        END DO
30      CONTINUE
C
101     FORMAT(A,$)
        END
C
C******************************************************************************
C Calcula el maximo de una linea alrededor del canal I, devolviendo la
C posicion del maximo calculada a traves del ajuste a una parabola.
        REAL FUNCTION FINDMAX(I)
        IMPLICIT NONE
        INTEGER I
C
        INCLUDE 'redlib.inc'
        INTEGER K
        INTEGER NSIDE
        INTEGER NFIT
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL SX(NCMAX),SY(NCMAX)
        REAL XF(NCMAX),YF(NCMAX)
        REAL A(3),CHISQR
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKNCHAN/NCHAN
        COMMON/BLKS/SX,SY
        COMMON/BLKNSIDE/NSIDE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO K=I-NSIDE,I+NSIDE
          XF(K-I+NSIDE+1)=SX(K)
          YF(K-I+NSIDE+1)=SY(K)
        END DO
        NFIT=2*NSIDE+1
C
        CALL POLFIT(XF,YF,YF,NFIT,3,0,A,CHISQR)
C
        DO K=1,101
          XF(K)=SX(I-NSIDE)+REAL(2*NSIDE)*REAL(K-1)/100
          YF(K)=A(1)+A(2)*XF(K)+A(3)*XF(K)*XF(K)
        END DO
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGLINE(101,XF,YF)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
C
        FINDMAX=-A(2)/(2.0*A(3))
C
        END
