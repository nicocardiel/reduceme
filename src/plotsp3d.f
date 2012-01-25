C------------------------------------------------------------------------------
C Version 18-September-1997                                      file: plotsp3d
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
C Program: plotsp3d
C Classification: graphic display
C Description: Plots successive spectra with a 3-D perspective.
C
Comment
C
C Dibuja espectros en profundidad para dar sensacion de imagen tridimensional
C
        PROGRAM PLOTSP3D
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        REAL READF
C
        INTEGER NBUT
        PARAMETER(NBUT=6)
C
        INTEGER NB
        INTEGER I,J,KK,L
        INTEGER I0
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NNC1,NNC2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XOFFSET,YOFFSET
        REAL XMAX,XMIN,YMAX,YMAX0,YMIN,YMIN0
        REAL DX,DY
        REAL XC,YC
        REAL S(NCMAX),A(NCMAX,NSMAX)
        REAL X(NCMAX),F(NCMAX)
        REAL FMIN,FMAX
        REAL XV1,XV2,YV1,YV2
        CHARACTER*1 CH,COPC,CCONT
        CHARACTER*10 BLABEL(NBUT)
        CHARACTER*40 POSICION,CDUMMY
        CHARACTER*75 INFILE
        LOGICAL REPLOT
        LOGICAL LCOLOR(MAX_ID_RED)
C------------------------------------------------------------------------------
        DATA (BLABEL(KK),KK=1,NBUT)/'Up','Down','Zoom','Whole',
     +   'Y-cuts','EXIT'/
C
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='plotsp3d'
        CALL WELCOME('18-September-1997')
C
        CCONT='0'
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL BUTTQPR(XV1,XV2,YV1,YV2)
        CALL BUTTSPR(XV1,XV2,YV1,0.80)
        CALL BUTTQBR(XV1,XV2,YV1,YV2)
        CALL BUTTSBR(XV1,XV2,0.83,YV2)
        DO KK=1,NBUT
          CALL BUTTON(KK,BLABEL(KK),0)
        END DO
C
10      WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
20      WRITE(*,100)'Plot whole image (y/n) '
        COPC(1:1)=READC('y','yn')
        NS1=1
        NS2=NSCAN
        NC1=1
        NC2=NCHAN
        IF(COPC.EQ.'n')THEN
31        WRITE(CDUMMY,'(I10,A1,I10)')NC1,',',NC2
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'1st & last channel '
          CALL READ2I(CDUMMY(1:L),NC1,NC2)
          IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
            WRITE(*,101)'ERROR: invalid entry. Try again.'
            GOTO 31
          END IF
          WRITE(CDUMMY,'(I10,A1,I10)')NS1,',',NS2
          CALL RMBLANK(CDUMMY,CDUMMY,L)
32        WRITE(*,100)'1st & last scan    '
          CALL READ2I(CDUMMY(1:L),NS1,NS2)
          IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
            WRITE(*,101)'ERROR: invalid entry. Try again.'
            GOTO 32
          END IF
        END IF
C
        YMIN0=A(NC1,NS1)
        YMAX0=YMIN0
        DO I=NS1,NS2
          DO J=NC1,NC2
            YMIN0=AMIN1(YMIN0,A(J,I))
            YMAX0=AMAX1(YMAX0,A(J,I))
          END DO
        END DO
        WRITE(*,*)
        WRITE(*,100)'YMIN = '
        WRITE(*,*)YMIN0
        WRITE(*,100)'YMAX = '
        WRITE(*,*)YMAX0
        WRITE(*,100)'Change Limits (y/n) '
        COPC(1:1)=READC('n','yn')
        IF(COPC.EQ.'y')THEN
          WRITE(*,100)'YMIN '
          WRITE(CDUMMY,*)YMIN0
          YMIN0=READF(CDUMMY)
          WRITE(*,100)'YMAX '
          WRITE(CDUMMY,*)YMAX0
          YMAX0=READF(CDUMMY)
        END IF
C
        I0=(NS1+NS2)/2
C
35      CONTINUE
C        WRITE(*,100)'X-offset between each spectrum (channels)'
C        XOFFSET=READF('@')
        XOFFSET=REAL(NC2-NC1+1)/(4*REAL(NS2-NS1+1))
C        WRITE(*,100)'Y-offset between each spectrum (counts)  '
C        YOFFSET=READF('@')
        YOFFSET=0.8*(YMAX0-YMIN0)/REAL(NS2-NS1+1)
        XMIN=REAL(NC1)-1.
        XMAX=REAL(NC2)+REAL(NS2-NS1)*XOFFSET+1.
        YMAX=YMAX0+REAL(NS2-NS1-1)*YOFFSET
        YMAX=YMAX+YOFFSET
        YMIN=YMIN0-YOFFSET
        DY=YMAX-YMIN
        DX=XMAX-XMIN
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.83)
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          END IF
          CALL PGIDEN_RED
          CALL PGLABEL('channel','counts',CHAR(32))
          IF(TRUELEN(OBJECT).GT.0)THEN
            CALL PGPTEXT(XMIN,YMAX+DY/80,0.,0.,
     +       'file: '//INFILE(1:TRUELEN(INFILE))//
     +       ' ['//OBJECT(1:TRUELEN(OBJECT))//']')
          ELSE
            CALL PGPTEXT(XMIN,YMAX+DY/80,0.,0.,'file: '//INFILE)
          END IF
          IF(LCOLOR(ITERM))THEN
            CALL PGSCI(2)
          ELSE
            CALL PGSLS(2)
          END IF
        END DO
        DO I=NS2,NS1,-1
          DO J=NC1,NC2
            X(J-NC1+1)=REAL(J)+REAL(I-NS1)*XOFFSET
            S(J-NC1+1)=A(J,I)+REAL(I-NS1)*YOFFSET
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
          END DO
        END DO
        CH=' '
        DO J=NC1,NC2
          X(J-NC1+1)=REAL(J)+REAL(I0-NS1)*XOFFSET
          S(J-NC1+1)=A(J,I0)+REAL(I0-NS1)*YOFFSET
          F(J-NC1+1)=A(J,I0)
        END DO
        CALL FINDMM(NC2-NC1+1,F,FMIN,FMAX)
        CALL PRINTD(I0,FMIN,FMAX)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM))THEN
            CALL PGSCI(1)
          ELSE
            CALL PGSLS(1)
          END IF
          CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
          WRITE(POSICION,'(A,I4)')'Scan # ',I0
          CALL PGPTEXT(XMIN+DX/40,YMAX-DY/20,0.,0.,POSICION)
        END DO
C------------------------------------------------------------------------------
        DO WHILE(CH.NE.'Q')
          CH=' '
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
c..............................................................................
          IF(NB.EQ.0)THEN
            WRITE(*,100)'Channel: '
            WRITE(*,*)XC-REAL(I0-NS1)*XOFFSET
          ELSE
            CALL BUTTON(NB,BLABEL(NB),4)
            CALL BUTTON(NB,BLABEL(NB),1)
          END IF
c..............................................................................
          IF((NB.EQ.1).OR.(NB.EQ.2))THEN
            REPLOT=.TRUE.
          ELSE
            REPLOT=.FALSE.
          END IF
          IF(NB.EQ.1)THEN
            I0=I0+1
          END IF
          IF(NB.EQ.2)THEN
             I0=I0-1
          END IF
          IF(NB.EQ.6)THEN
            CH='Q'
          END IF
          IF(REPLOT)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM))THEN
                CALL PGSCI(2)
              ELSE
                CALL PGSCI(0)
              END IF
              CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
              CALL PGSCI(1)
              IF(.NOT.LCOLOR(ITERM))THEN
                CALL PGSLS(2)
                CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
                CALL PGSLS(1)
              END IF
            END DO
            IF(I0.LT.NS1) I0=NS2
            IF(I0.GT.NS2) I0=NS1
            DO J=NC1,NC2
              X(J-NC1+1)=REAL(J)+REAL(I0-NS1)*XOFFSET
              S(J-NC1+1)=A(J,I0)+REAL(I0-NS1)*YOFFSET
              F(J-NC1+1)=A(J,I0)
            END DO
            CALL FINDMM(NC2-NC1+1,F,FMIN,FMAX)
            CALL PRINTD(I0,FMIN,FMAX)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
              CALL PGSCI(0)
              CALL PGPTEXT(XMIN+DX/40,YMAX-DY/20,0.,0.,POSICION)
              CALL PGSCI(1)
              WRITE(POSICION,'(A,I4)')'Scan # ',I0
              CALL PGPTEXT(XMIN+DX/40,YMAX-DY/20,0.,0.,POSICION)
            END DO
          END IF
c..............................................................................
          IF(NB.EQ.3)THEN
            WRITE(*,101)'Press region to be zoomed (first and last '//
     +       'channel)...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            NNC1=INT(XC-REAL(I0-NS1)*XOFFSET+0.5)
            IF(NNC1.LT.1) NNC1=1
            IF(NNC1.GT.NCHAN) NNC1=NCHAN
            WRITE(*,110)'First channel: ',NNC1
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(4,0,XC,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            NNC2=INT(XC-REAL(I0-NS1)*XOFFSET+0.5)
            IF(NNC2.LT.1) NNC2=1
            IF(NNC2.GT.NCHAN) NNC2=NCHAN
            WRITE(*,110)'Last  channel: ',NNC2
            IF(NNC1.EQ.NNC2)THEN
              WRITE(*,101)'ERROR: Invalid limits.'
            ELSE
              NC1=NNC1
              NC2=NNC2
ccc              CALL RPGERASW(0.,1.,0.,0.83)
              CALL BUTTON(NB,BLABEL(NB),0)
              GOTO 35
            END IF
          END IF
c..............................................................................
          IF(NB.EQ.4)THEN
            NC1=1
            NC2=NCHAN
ccc            CALL RPGERASW(0.,1.,0.,0.83)
            CALL BUTTON(NB,BLABEL(NB),0)
            GOTO 35
          END IF
c..............................................................................
          IF(NB.EQ.5)THEN
            WRITE(*,100)'New YMIN '
            WRITE(CDUMMY,*)YMIN
            YMIN0=READF(CDUMMY)
            WRITE(*,100)'New YMAX '
            WRITE(CDUMMY,*)YMAX0
            YMAX0=READF(CDUMMY)
ccc            CALL RPGERASW(0.,1.,0.,0.83)
            CALL BUTTON(NB,BLABEL(NB),0)
            GOTO 35
          END IF
c..............................................................................
          IF(NB.NE.0) CALL BUTTON(NB,BLABEL(NB),0)
c..............................................................................
        END DO
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,101)'(1) Same image'
        WRITE(*,101)'(2) New image'
        WRITE(*,101)'(0) QUIT'
        WRITE(*,100)'Option '
        CCONT(1:1)=READC(CCONT,'012')
        IF(CCONT.EQ.'1')THEN
          GOTO 20
        ELSEIF(CCONT.EQ.'2')THEN
          GOTO 10
        ELSEIF(CCONT.EQ.'0')THEN
          CALL PGEND
        END IF
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE PRINTD(I0,FMIN,FMAX)
        IMPLICIT NONE
        INTEGER I0
        INTEGER L,M,N
        REAL FMIN,FMAX
        CHARACTER*70 CDUMMY,CDUMMYF
C
        WRITE(CDUMMY,*)'Scan#',I0
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CDUMMYF(1:L)=CDUMMY(1:L)
        WRITE(CDUMMY,*)'Min=',FMIN
        CALL RMBLANK(CDUMMY,CDUMMY,M)
        CDUMMYF(L+1:L+5)='     '
        CDUMMYF(L+1+5:L+M+5)=CDUMMY(1:M)
        WRITE(CDUMMY,*)'Max=',FMAX
        CALL RMBLANK(CDUMMY,CDUMMY,N)
        CDUMMYF(L+M+5+1:L+M+5+5)='     '
        CDUMMYF(L+M+5+5+1:L+M+5+5+N)=CDUMMY(1:N)
        WRITE(*,'(A)')CDUMMYF(1:L+M+N+10)
        END
