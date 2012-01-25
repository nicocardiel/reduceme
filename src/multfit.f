C------------------------------------------------------------------------------
C Version 2-May-1999                                             file:multfit.f
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
C Program: multfit
C Classification: arithmetic & manipulations
C Description: perform simultaneous fits.
C
Comment
C
C NOTA: para incluir nuevas funciones, buscar la cadena CADD y seguir
C instrucciones.
C
        PROGRAM MULTFIT
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXPARAM                   !numero maximo de parametros/perfil
        PARAMETER (NMAXPARAM=7)
        INTEGER NMAXFUNCT                            !numero maximo de perfiles
        PARAMETER (NMAXFUNCT=12)
        INTEGER NMAXT
        PARAMETER (NMAXT=NMAXPARAM*NMAXFUNCT)
C
        INTEGER I,J,L
        INTEGER NB,NBLOCAL
        INTEGER N1,N2,NTOT,NPT
        INTEGER IXC1,IXC2
        INTEGER NCMAX_LOCAL,NSMAX_LOCAL
        INTEGER NCBUFF,NBUFF
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NF1,NF2
        INTEGER NDIM,NEVAL
        INTEGER NPARAM(NMAXT)
        INTEGER NCOLOR
        REAL A(NCMAX,NSMAX)
        REAL X(NCMAX),Y(NCMAX),YPLOT(NCMAX),YPLOTT(NCMAX)
        REAL XMIN,XMAX,YMIN,YMAX,DX,DY
        REAL XC,YC
        REAL PARAM(NMAXPARAM,NMAXFUNCT)
        REAL EEPARAM(NMAXPARAM,NMAXFUNCT)
        REAL XX0(NMAXT),DXX0(NMAXT),XXF(NMAXT),DXXF(NMAXT)
        EXTERNAL FUNKMFIT
        REAL FUNKMFIT
        REAL YRMSTOL
        REAL POL
        REAL FCTEDOUBLET
        CHARACTER*1 CH,CSURE,CACTION
        CHARACTER*1 FIX(NMAXPARAM,NMAXFUNCT),CPROF(NMAXFUNCT)
        CHARACTER*1 FIX_TEMPORAL,CPROF_TEMPORAL
        CHARACTER*3 CFITB(NMAXFUNCT)
        CHARACTER*50 CDUMMY
        CHARACTER*75 FILENAME,XLABEL,YLABEL
        LOGICAL IFSCAN(NSMAX),IFCHAN(NCMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LDEFBUFF(NMAXFUNCT),LUSEBUFF(NMAXFUNCT)
        LOGICAL LFIRSTP,LDEFAULT,LANY,LREGION
        LOGICAL LBEXIST
C
        COMMON/BLKFITMF1/NPT,NF1,NF2
        COMMON/BLKFITMF2/X,Y
        COMMON/BLKFITMF3/NPARAM
        COMMON/BLKFITMF4/PARAM
        COMMON/BLKFITMF5/FIX,CPROF
        COMMON/BLKFITMF6/LUSEBUFF
        COMMON/BLKFITMF7/FCTEDOUBLET
C------------------------------------------------------------------------------
        DATA(CFITB(I),I=1,NMAXFUNCT)/'#1','#2','#3','#4','#5',
     +   '#6','#7','#8','#9','#10','#11','#12'/
C
        OUTFILEX=OUTFILEX
C 
        LFIRSTP=.TRUE.
        YLABEL='No. of counts'
        YRMSTOL=1.E-6
        LREGION=.FALSE.
C
        NCBUFF=0
        DO NBUFF=1,NMAXFUNCT
          LDEFBUFF(NBUFF)=.FALSE.
          LUSEBUFF(NBUFF)=.FALSE.
          DO I=1,NMAXPARAM
            PARAM(I,NBUFF)=0.
            FIX(I,NBUFF)='n'
          END DO
        END DO
C------------------------------------------------------------------------------
        THISPROGRAM='multfit'
        CALL WELCOME('3-May-1999')
C------------------------------------------------------------------------------
        NCMAX_LOCAL=NCMAX
        NSMAX_LOCAL=NSMAX
        IF(NCMAX_LOCAL.LT.NSMAX_LOCAL)THEN
          WRITE(*,101)'FATAL ERROR: this program requires '//
     +     'NCMAX.GT.NSMAX'
          STOP
        END IF
C------------------------------------------------------------------------------
C abrimos salida grafica
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL BUTTSYB(4)
C------------------------------------------------------------------------------
        WRITE(*,100)'Input file name'
        FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        WRITE(*,100)'Reading file...'
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        WRITE(*,101)'  ...OK! File read and closed.'
C------------------------------------------------------------------------------
        CALL BUTTON(1,'[x] cut',0)
        CALL BUTTON(2,'[y] cut',0)
        CALL BUTTON(3,'[z]oom',0)
        CALL BUTTON(3,'[z]oom',3)
        CALL BUTTON(4,'[w]hole',0)
        CALL BUTTON(4,'[w]hole',3)
        CALL BUTTON(5,'Y[-]limits',0)
        CALL BUTTON(5,'Y[-]limits',3)
        CALL BUTTON(6,'[q]uit',0)
        CALL BUTTON(7,'[p]rofile',0)
        CALL BUTTON(7,'[p]rofile',3)
        CALL BUTTON(8,'yrms[t]ol',0)
        CALL BUTTON(8,'yrms[t]ol',3)
        CALL BUTTON(9,'[r]egion',0)
        CALL BUTTON(9,'[r]egion',3)
        CALL BUTTON(10,'[m]fit',0)
        CALL BUTTON(10,'[m]fit',3)
        CALL BUTTON(11,'[l]ist',0)
        CALL BUTTON(11,'[l]ist',3)
        CALL BUTTON(12,'r[e]plot',0)
        CALL BUTTON(12,'r[e]plot',3)
        DO I=1,12
          CALL BUTTON(I+12,CFITB(I),0)
          CALL BUTTON(I+12,CFITB(I),3)
        END DO
C------------------------------------------------------------------------------
7       CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        NBLOCAL=INDEX('xyzw-qptrmle',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.0)THEN
          WRITE(*,100)'Cursor at (x,y): '
          WRITE(CDUMMY,*) XC
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'('
          WRITE(*,100) CDUMMY(1:L)
          WRITE(*,100)','
          WRITE(CDUMMY,*) YC
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100) CDUMMY(1:L)
          WRITE(*,101)')'
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.1)THEN
          CALL BUTTON(1,'[x] cut',1)
          CALL BUTTON(2,'[y] cut',0)
          DO I=1,NSCAN
            IFSCAN(I)=.FALSE.
          END DO
10        WRITE(*,100)'Scan region (0,0=EXIT) '
          CALL READ2I('0,0',N1,N2)
          IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 12
          IF((N1.LT.1).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
            GOTO 10
          END IF
          DO I=N1,N2
            IFSCAN(I)=.TRUE.
          END DO
          GOTO 10
12        NTOT=0
          DO I=1,NSCAN
            IF(IFSCAN(I)) NTOT=NTOT+1
          END DO
          IF(NTOT.EQ.0)THEN
            WRITE(*,101)'ERROR: number of scans added = 0. Try again.'
            GOTO 10
          ELSE
            WRITE(*,100)'>>> No. of scans added: '
            WRITE(*,*) NTOT
          END IF
          DO J=1,NCHAN
            X(J)=REAL(J)
            Y(J)=0.
            DO I=1,NSCAN
              IF(IFSCAN(I)) Y(J)=Y(J)+A(J,I)
            END DO
            Y(J)=Y(J)/REAL(NTOT)
          END DO
          XLABEL='channel'
          NPT=NCHAN
          N1=1
          N2=NPT
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'[y] cut',1)
          CALL BUTTON(1,'[x] cut',0)
          DO J=1,NCHAN
            IFCHAN(J)=.FALSE.
          END DO
20        WRITE(*,100)'Channel region (0,0=EXIT) '
          CALL READ2I('0,0',N1,N2)
          IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 22
          IF((N1.LT.1).OR.(N2.GT.NCHAN).OR.(N1.GT.N2))THEN
            WRITE(*,101)'ERROR: numbers out of range. Try again.'
            GOTO 20
          END IF
          DO J=N1,N2
            IFCHAN(J)=.TRUE.
          END DO
          GOTO 20
22        NTOT=0
          DO J=1,NCHAN
            IF(IFCHAN(J)) NTOT=NTOT+1
          END DO
          IF(NTOT.EQ.0)THEN
            WRITE(*,101)'ERROR: number of channels added = 0. '//
     +       'Try again.'
            GOTO 20
          ELSE
            WRITE(*,100)'>>> No. of channels added: '
            WRITE(*,*) NTOT
          END IF
          DO I=1,NSCAN
            X(I)=REAL(I)
            Y(I)=0.
            DO J=1,NCHAN
              IF(IFCHAN(J)) Y(I)=Y(I)+A(J,I)
            END DO
            Y(I)=Y(I)/REAL(NTOT)
          END DO
          XLABEL='scan'
          NPT=NSCAN
          N1=1
          N2=NPT
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[z]oom',5)
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IXC1=INT(XC+0.5)
          IF(IXC1.LT.N1) IXC1=N1
          IF(IXC1.GT.N2) IXC1=N2
          WRITE(*,*)
          WRITE(*,100)'Point #1, cursor at x= '
          WRITE(*,*)INT(XC+0.5)
          WRITE(*,100)'Press mouse button...'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(4,0,REAL(IXC1),0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IXC2=INT(XC+0.5)
          IF(IXC2.LT.N1) IXC2=N1
          IF(IXC2.GT.N2) IXC2=N2
          WRITE(*,*)
          WRITE(*,100)'Point #2, cursor at x= '
          WRITE(*,*)INT(XC+0.5)
          N1=MIN0(IXC1,IXC2)
          N2=MAX0(IXC1,IXC2)
          CALL BUTTON(3,'[z]oom',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,'[w]hole',5)
          N1=1
          N2=NPT
          CALL BUTTON(4,'[w]hole',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.5)THEN
          CALL BUTTON(5,'Y[-]limits',5)
          WRITE(CDUMMY,*)YMIN
          WRITE(*,100)'Ymin '
          YMIN=READF(CDUMMY)
          WRITE(CDUMMY,*)YMAX
          WRITE(*,100)'Ymax '
          YMAX=READF(CDUMMY)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              IF(LFIRSTP)THEN
                LFIRSTP=.FALSE.
              ELSE
                CALL RPGERASW(0.,1.,0.,0.80)
              END IF
              CALL PGSLW(3)
              CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
              CALL PGSLW(1)
              CALL PGIDEN_RED
            ELSE
              CALL PGSLW(3)
              CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
              CALL PGSLW(1)
              CALL PGIDEN_RED
            END IF
            CALL PGLABEL(XLABEL,YLABEL,FILENAME)
            CALL PGBIN(NPT,X,Y,.TRUE.)
            IF(LREGION)THEN
              IF(LCOLOR(ITERM)) CALL PGSCI(7)
              CALL PGSLS(4)
              CALL PGMOVE(REAL(NF1),YMIN)
              CALL PGDRAW(REAL(NF1),YMAX)
              CALL PGMOVE(REAL(NF2),YMIN)
              CALL PGDRAW(REAL(NF2),YMAX)
              CALL PGSLS(1)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            ELSE
              NF1=N1
              NF2=N2
            END IF
          END DO
          CALL BUTTON(5,'Y[-]limits',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[q]uit',5)
          WRITE(*,100)'Do you really want to quit (y/n) '
          CSURE(1:1)=READC('y','yn')
          IF(CSURE.EQ.'y') GOTO 90
          CALL BUTTON(6,'[q]uit',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,'[p]rofile',5)
          CALL BUTTON(8,'yrms[t]ol',0)
          CALL BUTTON(9,'[r]egion',0)
          CALL BUTTON(11,'[l]ist',0)
          CALL BUTTON(12,'r[e]plot',0)
C
          WRITE(*,150)
          WRITE(*,100)'Buffer number '
          NCBUFF=READILIM('@',1,NMAXFUNCT)
          IF(.NOT.LDEFBUFF(NCBUFF))THEN
            LDEFAULT=.FALSE.
            CALL BUTTON(NCBUFF+12,CFITB(NCBUFF),1)
            CALL BUTTON(10,'[m]fit',2)
            LDEFBUFF(NCBUFF)=.TRUE.
            LUSEBUFF(NCBUFF)=.TRUE.
CADD: introducir funcion en menu 
            WRITE(*,101)'-----------------------------------'
            WRITE(*,101)'(c) Cauchy function'
            WRITE(*,101)'(g) Gaussian'
            WRITE(*,101)'(d) two related gaussians'
            WRITE(*,101)'(p) Polynomial'
            WRITE(*,101)'-----------------------------------'
            WRITE(*,100)'Option '
            CPROF_TEMPORAL(1:1)=READC('@','cgpd')
            CPROF(NCBUFF)=CPROF_TEMPORAL
          ELSE
            WRITE(*,101)'WARNING: This buffer has already been defined.'
            WRITE(*,100)'[m]odify or [d]elete this buffer (d/m) '
            CACTION(1:1)=READC('m','md')
            IF(CACTION.EQ.'d')THEN
              LDEFBUFF(NCBUFF)=.FALSE.
              LUSEBUFF(NCBUFF)=.FALSE.
              CALL BUTTON(NCBUFF+12,CFITB(NCBUFF),3)
              LANY=.FALSE.
              DO I=1,NMAXFUNCT
                IF(LDEFBUFF(I)) LANY=.TRUE.
              END DO
              IF(.NOT.LANY)THEN
                CALL BUTTON(8,'yrms[t]ol',3)
                CALL BUTTON(9,'[r]egion',3)
                CALL BUTTON(10,'[m]fit',3)
                CALL BUTTON(11,'[l]ist',3)
                CALL BUTTON(12,'r[e]plot',3)
              END IF
              CALL BUTTON(7,'[p]rofile',0)
              GOTO 7
            END IF
            LDEFAULT=.TRUE.
          END IF
          WRITE(*,*)
C
CADD: incluir nueva entrada en estructura IF
c..............................................................................
          IF(CPROF(NCBUFF).EQ.'c')THEN
            WRITE(*,101)'Cauchy function: y=amp/[sigma^2+(x-x0)^2]'
            NPARAM(NCBUFF)=3
c
            WRITE(*,100)'> amp.......'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(1,NCBUFF)
              PARAM(1,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(1,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(1,NCBUFF),'yn')
              FIX(1,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(1,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(1,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(1,NCBUFF)
                EEPARAM(1,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(1,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> sigma.....'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(2,NCBUFF)
              PARAM(2,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(2,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(2,NCBUFF),'yn')
              FIX(2,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(2,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(2,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(2,NCBUFF)
                EEPARAM(2,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(2,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> x0........'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(3,NCBUFF)
              PARAM(3,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(3,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(3,NCBUFF),'yn')
              FIX(3,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(3,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(3,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(3,NCBUFF)
                EEPARAM(3,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(3,NCBUFF)=READF('@')
              END IF
            END IF
c..............................................................................
          ELSEIF(CPROF(NCBUFF).EQ.'g')THEN
            WRITE(*,101)'Gaussian: y=amp*exp[-(x-x0)^2/(2*sigma^2)]'
            NPARAM(NCBUFF)=3
c
            WRITE(*,100)'> amp.......'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(1,NCBUFF)
              PARAM(1,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(1,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(1,NCBUFF),'yn')
              FIX(1,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(1,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(1,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(1,NCBUFF)
                EEPARAM(1,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(1,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> sigma.....'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(2,NCBUFF)
              PARAM(2,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(2,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(2,NCBUFF),'yn')
              FIX(2,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(2,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(2,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(2,NCBUFF)
                EEPARAM(2,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(2,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> x0........'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(3,NCBUFF)
              PARAM(3,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(3,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(3,NCBUFF),'yn')
              FIX(3,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(3,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(3,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(3,NCBUFF)
                EEPARAM(3,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(3,NCBUFF)=READF('@')
              END IF
            END IF
c..............................................................................
          ELSEIF(CPROF(NCBUFF).EQ.'d')THEN
            WRITE(*,101)'Gaussian: y=amp1*exp[-(x-x01)^2/(2*sigma1^2)]'
            WRITE(*,101)'           +amp2*exp[-(x-x02)^2/(2*sigma2^2)]'
            WRITE(*,100)'where area2=cte*area1 <==> '
            WRITE(*,101)'amp2=amp1*sigma1/sigma2*cte'
            NPARAM(NCBUFF)=5
c
            WRITE(*,100)'> Constant relation between areas '
            FCTEDOUBLET=READF('3.0')
c
            WRITE(*,100)'> amp1......'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(1,NCBUFF)
              PARAM(1,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(1,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(1,NCBUFF),'yn')
              FIX(1,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(1,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(1,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(1,NCBUFF)
                EEPARAM(1,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(1,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> sigma1....'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(2,NCBUFF)
              PARAM(2,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(2,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(2,NCBUFF),'yn')
              FIX(2,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(2,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(2,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(2,NCBUFF)
                EEPARAM(2,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(2,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> x01.......'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(3,NCBUFF)
              PARAM(3,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(3,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(3,NCBUFF),'yn')
              FIX(3,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(3,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(3,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(3,NCBUFF)
                EEPARAM(3,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(3,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> sigma2....'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(4,NCBUFF)
              PARAM(4,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(4,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(4,NCBUFF),'yn')
              FIX(4,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(4,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(4,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(4,NCBUFF)
                EEPARAM(4,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(4,NCBUFF)=READF('@')
              END IF
            END IF
c
            WRITE(*,100)'> x02.......'
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)PARAM(5,NCBUFF)
              PARAM(5,NCBUFF)=READF(CDUMMY)
            ELSE
              PARAM(5,NCBUFF)=READF('@')
            END IF
            WRITE(*,100)'fix this parameter (y/n) '
            IF(LDEFAULT)THEN
              FIX_TEMPORAL(1:1)=READC(FIX(5,NCBUFF),'yn')
              FIX(5,NCBUFF)=FIX_TEMPORAL
            ELSE
              FIX_TEMPORAL(1:1)=READC('n','yn')
              FIX(5,NCBUFF)=FIX_TEMPORAL
            END IF
            IF(FIX(5,NCBUFF).EQ.'n')THEN
              WRITE(*,100)'> length-scale '
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)EEPARAM(5,NCBUFF)
                EEPARAM(5,NCBUFF)=READF(CDUMMY)
              ELSE
                EEPARAM(5,NCBUFF)=READF('@')
              END IF
            END IF
c..............................................................................
          ELSEIF(CPROF(NCBUFF).EQ.'p')THEN
            WRITE(*,101)'Polynomial: y = a0+a1*x+a2*x*x+...'
            WRITE(*,100)'polynomial degree '
            IF(LDEFAULT)THEN
              WRITE(CDUMMY,*)NPARAM(NCBUFF)-1
              NPARAM(NCBUFF)=READILIM(CDUMMY,0,6)
              NPARAM(NCBUFF)=NPARAM(NCBUFF)+1
            ELSE
              NPARAM(NCBUFF)=READILIM('@',0,6)
              NPARAM(NCBUFF)=NPARAM(NCBUFF)+1
            END IF
            DO I=1,NPARAM(NCBUFF)
              WRITE(*,'(A3,I1,A8,$)')'> a',I-1,'........'
              IF(LDEFAULT)THEN
                WRITE(CDUMMY,*)PARAM(I,NCBUFF)
                PARAM(I,NCBUFF)=READF(CDUMMY)
              ELSE
                PARAM(I,NCBUFF)=READF('@')
              END IF
              WRITE(*,100)'fix this parameter (y/n) '
              IF(LDEFAULT)THEN
                FIX_TEMPORAL(1:1)=READC(FIX(I,NCBUFF),'yn')
                FIX(I,NCBUFF)=FIX_TEMPORAL
              ELSE
                FIX_TEMPORAL(1:1)=READC('n','yn')
                FIX(I,NCBUFF)=FIX_TEMPORAL
              END IF
              IF(FIX(I,NCBUFF).EQ.'n')THEN
                WRITE(*,100)'> length-scale '
                IF(LDEFAULT)THEN
                  WRITE(CDUMMY,*)EEPARAM(I,NCBUFF)
                  EEPARAM(I,NCBUFF)=READF(CDUMMY)
                ELSE
                  EEPARAM(I,NCBUFF)=READF('@')
                END IF
              END IF
            END DO
c..............................................................................
          END IF
          CALL BUTTON(7,'[p]rofile',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,'yrms[t]ol',5)
          WRITE(CDUMMY,*)YRMSTOL
          WRITE(*,100)'YRMSTOL for DOWNHILL '
          YRMSTOL=READF(CDUMMY)
          CALL BUTTON(8,'yrms[t]ol',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.9)THEN
          WRITE(*,101)'* Define region to be fitted: '
          WRITE(*,100)'N1 '
          WRITE(CDUMMY,*) NF1
          NF1=READILIM(CDUMMY,1,NPT)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(7)
            CALL PGSLS(4)
            CALL PGMOVE(REAL(NF1),YMIN)
            CALL PGDRAW(REAL(NF1),YMAX)
            CALL PGSLS(1)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          WRITE(*,100)'N2 '
          WRITE(CDUMMY,*) NF2
          NF2=READILIM(CDUMMY,NF1,NPT)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(7)
            CALL PGSLS(4)
            CALL PGMOVE(REAL(NF2),YMIN)
            CALL PGDRAW(REAL(NF2),YMAX)
            CALL PGSLS(1)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          LREGION=.TRUE.
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.10)THEN
          CALL BUTTON(10,'[m]fit',5)
          IF(.NOT.LREGION)THEN
            WRITE(*,101)'ERROR: define region to be employed first!'
            WRITE(*,100)'(press <CR> to continue...)'
            READ(*,*)
            CALL BUTTON(10,'[m]fit',0)
            GOTO 7
          END IF
C
          NDIM=0
          DO I=1,NMAXFUNCT
            IF(LUSEBUFF(I))THEN
              DO J=1,NPARAM(I)
                IF(FIX(J,I).EQ.'n')THEN
                  NDIM=NDIM+1
                  XX0(NDIM)=PARAM(J,I)
                  DXX0(NDIM)=EEPARAM(J,I)
                END IF
              END DO
            END IF
          END DO
          IF(NDIM.GT.0)THEN
            WRITE(*,150)
            WRITE(*,101)'* Parameters before DOWNHILL:'
            DO J=1,NDIM
              WRITE(*,*)J,XX0(J),DXX0(J)
            END DO
            DO J=1,NDIM
              IF(DXX0(J).EQ.0.0)THEN
                WRITE(*,101)'ERROR: length-scale value.EQ.0.0'
                WRITE(*,100)'> Paramer #'
                WRITE(*,*)J
                WRITE(*,100)'> Value..........:'
                WRITE(*,*)XX0(J)
                WRITE(*,100)'> length-scale...:'
                WRITE(*,*)DXX0(J)
                WRITE(*,100)'New length-scale value'
                DXX0(J)=READF('@')
              END IF
            END DO
            CALL DOWNHILL(NDIM,XX0,DXX0,FUNKMFIT,1.0,0.5,2.0,YRMSTOL,
     +       XXF,DXXF,NEVAL)
            WRITE(*,150)
            WRITE(*,101)'* Parameters after DOWNHILL (+rms):'
            DO J=1,NDIM
              WRITE(*,*)J,XXF(J),DXXF(J)
            END DO
            WRITE(*,150)
          ELSE
            WRITE(*,150)
            WRITE(*,101)'ERROR: number of parameters to be fitted = 0!'
            WRITE(*,100)'(press <CR> to continue...)'
            READ(*,*)
            WRITE(*,150)
          END IF
C
          NDIM=0
          DO J=N1,N2
            YPLOTT(J)=0.
          END DO
          NCOLOR=1
          DO I=1,NMAXFUNCT
            IF(LUSEBUFF(I))THEN
              DO J=1,NPARAM(I)
                IF(FIX(J,I).EQ.'n')THEN
                  NDIM=NDIM+1
                  PARAM(J,I)=XXF(NDIM)
                  EEPARAM(J,I)=DXXF(NDIM)
                END IF
              END DO
CADD: introducir funcion para dibujar el ajuste calculado por DOWNHILL
c..............................................................................
              IF(CPROF(I).EQ.'c')THEN
                DO J=N1,N2
                  YPLOT(J)=PARAM(1,I)/(PARAM(2,I)*PARAM(2,I)+
     >             (X(J)-PARAM(3,I))*(X(J)-PARAM(3,I)))
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              ELSEIF(CPROF(I).EQ.'g')THEN
                DO J=N1,N2
                  YPLOT(J)=PARAM(1,I)*
     >             EXP(-(X(J)-PARAM(3,I))*(X(J)-PARAM(3,I))/
     >             (2.*PARAM(2,I)*PARAM(2,I)))
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              ELSEIF(CPROF(I).EQ.'d')THEN
                DO J=N1,N2
                  YPLOT(J)=PARAM(1,I)*
     >             EXP(-(X(J)-PARAM(3,I))*(X(J)-PARAM(3,I))/
     >             (2.*PARAM(2,I)*PARAM(2,I)))
                  YPLOT(J)=YPLOT(J)+
     >             FCTEDOUBLET*PARAM(1,I)*PARAM(2,I)/PARAM(4,I)*
     >             EXP(-(X(J)-PARAM(5,I))*(X(J)-PARAM(5,I))/
     >             (2.*PARAM(4,I)*PARAM(4,I)))
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              ELSEIF(CPROF(I).EQ.'p')THEN
                DO J=N1,N2
                  POL=PARAM(NPARAM(I),I)
                  DO L=NPARAM(I)-1,1,-1
                    POL=POL*X(J)+PARAM(L,I)
                  END DO
                  YPLOT(J)=POL
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                NCOLOR=NCOLOR+1
                IF(NCOLOR.GT.14) NCOLOR=2
                IF(NCOLOR.EQ.3) NCOLOR=4
                IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
                IF(ITERM.EQ.1) CALL BUTTON(I+12,CFITB(I),-NCOLOR-1)
                CALL PGSLS(4)
                CALL PGBIN(NPT,X,YPLOT,.TRUE.)
                CALL PGSLS(1)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
            END IF
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGSLS(2)
            CALL PGBIN(NPT,X,YPLOTT,.TRUE.)
            CALL PGSLS(1)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
C
          CALL BUTTON(10,'[m]fit',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(11,'[l]ist',5)
          WRITE(*,150)
          WRITE(*,100)'Still undefined: '
          DO I=1,NMAXFUNCT
            IF(.NOT.LDEFBUFF(I))THEN
              IF(I.LT.10)THEN
                WRITE(*,'(A2,I1,$)')' #',I
              ELSE
                WRITE(*,'(A2,I2,$)')' #',I
              END IF
            END IF
          END DO
          WRITE(*,*)
          DO I=1,NMAXFUNCT
            IF(LDEFBUFF(I))THEN
              WRITE(*,'(A8,I2,A2,$)')'Buffer #',I,': '
CADD: incluir nombre de funcion y parametros en la estructura IF
              IF(CPROF(I).EQ.'c')THEN
                WRITE(*,101)'Cauchy function: y=amp/[sigma^2+(x-x0)^2]'
                WRITE(*,100)' > amp.......: '
                CALL MUESTRAPAR(FIX(1,I),PARAM(1,I))
                WRITE(*,100)' > sigma.....: '
                CALL MUESTRAPAR(FIX(2,I),PARAM(2,I))
                WRITE(*,100)' > x0........: '
                CALL MUESTRAPAR(FIX(3,I),PARAM(3,I))
              ELSEIF(CPROF(I).EQ.'g')THEN
                WRITE(*,101)'Gaussian: y=amp*exp[-(x-x0)^2/(2*sigma^2)]'
                WRITE(*,100)' > amp.......: '
                CALL MUESTRAPAR(FIX(1,I),PARAM(1,I))
                WRITE(*,100)' > sigma.....: '
                CALL MUESTRAPAR(FIX(2,I),PARAM(2,I))
                WRITE(*,100)' > x0........: '
                CALL MUESTRAPAR(FIX(3,I),PARAM(3,I))
              ELSEIF(CPROF(I).EQ.'d')THEN
                WRITE(*,100)'2 related gaussians: '
                WRITE(*,101)'y=amp1*exp[-(x-x01)^2/(2*sigma1^2)]'
                WRITE(*,100)'                     '
                WRITE(*,101)' +amp2*exp[-(x-x02)^2/(2*sigma2^2)]'
                WRITE(*,101)'where area2=cte*area1'
                WRITE(*,100)' > amp1......: '
                CALL MUESTRAPAR(FIX(1,I),PARAM(1,I))
                WRITE(*,100)' > sigma1....: '
                CALL MUESTRAPAR(FIX(2,I),PARAM(2,I))
                WRITE(*,100)' > x01.......: '
                CALL MUESTRAPAR(FIX(3,I),PARAM(3,I))
                WRITE(*,100)' > sigma2....: '
                CALL MUESTRAPAR(FIX(4,I),PARAM(4,I))
                WRITE(*,100)' > x02.......: '
                CALL MUESTRAPAR(FIX(5,I),PARAM(5,I))
              ELSEIF(CPROF(I).EQ.'p')THEN
                WRITE(*,101)'Polynomial: y = a0+a1*x+a2*x*x+...'
                WRITE(*,'(A,I1,A)')'(degree: ',NPARAM(I)-1,')'
                DO L=1,NPARAM(I)
                  WRITE(*,'(A4,I1,A,$)')' > a',L-1,'........: '
                  CALL MUESTRAPAR(FIX(L,I),PARAM(L,I))
                END DO
              END IF
            END IF
          END DO
          CALL BUTTON(11,'[l]ist',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(12,'r[e]plot',5)
          NDIM=0
          DO J=N1,N2
            YPLOTT(J)=0.
          END DO
          NCOLOR=1
          DO I=1,NMAXFUNCT
            IF(LUSEBUFF(I))THEN
CADD: introducir funcion para dibujar el ajuste calculado por DOWNHILL
C     (NOTA: el codigo aqui es identico que en NB.EQ.10  --ver mas arriba--)
c..............................................................................
              IF(CPROF(I).EQ.'c')THEN
                DO J=N1,N2
                  YPLOT(J)=PARAM(1,I)/(PARAM(2,I)*PARAM(2,I)+
     >             (X(J)-PARAM(3,I))*(X(J)-PARAM(3,I)))
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              ELSEIF(CPROF(I).EQ.'g')THEN
                DO J=N1,N2
                  YPLOT(J)=PARAM(1,I)*
     >             EXP(-(X(J)-PARAM(3,I))*(X(J)-PARAM(3,I))/
     >             (2.*PARAM(2,I)*PARAM(2,I)))
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              ELSEIF(CPROF(I).EQ.'d')THEN
                DO J=N1,N2
                  YPLOT(J)=PARAM(1,I)*
     >             EXP(-(X(J)-PARAM(3,I))*(X(J)-PARAM(3,I))/
     >             (2.*PARAM(2,I)*PARAM(2,I)))
                  YPLOT(J)=YPLOT(J)+
     >             FCTEDOUBLET*PARAM(1,I)*PARAM(2,I)/PARAM(4,I)*
     >             EXP(-(X(J)-PARAM(5,I))*(X(J)-PARAM(5,I))/
     >             (2.*PARAM(4,I)*PARAM(4,I)))
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              ELSEIF(CPROF(I).EQ.'p')THEN
                DO J=N1,N2
                  POL=PARAM(NPARAM(I),I)
                  DO L=NPARAM(I)-1,1,-1
                    POL=POL*X(J)+PARAM(L,I)
                  END DO
                  YPLOT(J)=POL
                  YPLOTT(J)=YPLOTT(J)+YPLOT(J)
                END DO
c..............................................................................
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                NCOLOR=NCOLOR+1
                IF(NCOLOR.GT.14) NCOLOR=2
                IF(NCOLOR.EQ.3) NCOLOR=4
                IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR)
                IF(ITERM.EQ.1) CALL BUTTON(I+12,CFITB(I),-NCOLOR-1)
                CALL PGSLS(4)
                CALL PGBIN(NPT,X,YPLOT,.TRUE.)
                CALL PGSLS(1)
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
            END IF
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGSLS(2)
            CALL PGBIN(NPT,X,YPLOTT,.TRUE.)
            CALL PGSLS(1)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
          CALL BUTTON(12,'r[e]plot',0)
C------------------------------------------------------------------------------
        ELSEIF((NB.GE.13).AND.(NB.LE.24))THEN
          IF(LUSEBUFF(NB-12))THEN
            LUSEBUFF(NB-12)=.FALSE.
            CALL BUTTON(NB,CFITB(NB-12),2)
          ELSE
            LUSEBUFF(NB-12)=.TRUE.
            CALL BUTTON(NB,CFITB(NB-12),1)
          END IF
          LANY=.FALSE.                          !vemos si queda algun buffer ON
          DO NBUFF=1,NMAXFUNCT
            IF(LUSEBUFF(NBUFF)) LANY=.TRUE.
          END DO
          IF(LANY)THEN
            CALL BUTTON(10,'[m]fit',2)
          ELSE
            CALL BUTTON(10,'[m]fit',3)
          END IF
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        IF((NB.GE.1).AND.(NB.LE.4))THEN
          CALL FINDMML(NPT,N1,N2,X,XMIN,XMAX)
          DX=XMAX-XMIN
          XMIN=XMIN-DX/20.
          XMAX=XMAX+DX/20.
          CALL BUTTON(3,'[z]oom',0)
          CALL BUTTON(4,'[w]hole',0)
          CALL BUTTON(5,'Y[-]limits',0)
          CALL BUTTON(7,'[p]rofile',0)
          CALL FINDMML(NPT,N1,N2,Y,YMIN,YMAX)
          DY=YMAX-YMIN
          YMIN=YMIN-DY/20.
          YMAX=YMAX+DY/20.
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              IF(LFIRSTP)THEN
                LFIRSTP=.FALSE.
              ELSE
                CALL RPGERASW(0.,1.,0.,0.80)
              END IF
              CALL PGSLW(3)
              CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
              CALL PGSLW(1)
              CALL PGIDEN_RED
            ELSE
              CALL PGSLW(3)
              CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
              CALL PGSLW(1)
              CALL PGIDEN_RED
            END IF
            CALL PGLABEL(XLABEL,YLABEL,FILENAME)
            CALL PGBIN(NPT,X,Y,.TRUE.)
            IF(LREGION)THEN
              IF(LCOLOR(ITERM)) CALL PGSCI(7)
              CALL PGSLS(3)
              CALL PGMOVE(REAL(NF1),YMIN)
              CALL PGDRAW(REAL(NF1),YMAX)
              CALL PGMOVE(REAL(NF2),YMIN)
              CALL PGDRAW(REAL(NF2),YMAX)
              CALL PGSLS(1)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            ELSE
              NF1=N1
              NF2=N2
            END IF
          END DO
        END IF
        GOTO 7
C------------------------------------------------------------------------------
90     CALL PGEND
        STOP
C
100     FORMAT(A,$)
101     FORMAT(A)
150     FORMAT(79('-'))
        END
C
C******************************************************************************
C
        SUBROUTINE MUESTRAPAR(FIX,PARAM)
        IMPLICIT NONE
        CHARACTER*1 FIX
        REAL PARAM
C
        INTEGER L
        CHARACTER*50 CDUMMY
C-----------------------------------------------------------------------------
        WRITE(CDUMMY,*)PARAM
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)CDUMMY(1:L)
        IF(FIX.EQ.'y')THEN
          WRITE(*,101)' <FIXED>'
        ELSE
          WRITE(*,101)' '
        END IF
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        REAL FUNCTION FUNKMFIT(XX)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER NMAXPARAM                   !numero maximo de parametros/perfil
        PARAMETER (NMAXPARAM=7)
        INTEGER NMAXFUNCT                            !numero maximo de perfiles
        PARAMETER (NMAXFUNCT=12)
        INTEGER NMAXT
        PARAMETER (NMAXT=NMAXPARAM*NMAXFUNCT)
C
        REAL XX(NMAXT)
C
        INTEGER I,J,K
        INTEGER NPT,NF1,NF2
        INTEGER NPARAM(NMAXT),NDIM
        REAL FF
        REAL X(NCMAX),Y(NCMAX)
        REAL YP(NCMAX,NMAXFUNCT)
        REAL PARAM(NMAXPARAM,NMAXFUNCT)
        REAL LPAR(NMAXPARAM)
        REAL FCTEDOUBLET
        CHARACTER*1 FIX(NMAXPARAM,NMAXFUNCT),CPROF(NMAXFUNCT)
        DOUBLE PRECISION DFUNK
        LOGICAL LUSEBUFF(NMAXFUNCT)
C
        COMMON/BLKFITMF1/NPT,NF1,NF2
        COMMON/BLKFITMF2/X,Y
        COMMON/BLKFITMF3/NPARAM
        COMMON/BLKFITMF4/PARAM
        COMMON/BLKFITMF5/FIX,CPROF
        COMMON/BLKFITMF6/LUSEBUFF
        COMMON/BLKFITMF7/FCTEDOUBLET
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        NDIM=0
        DO I=1,NMAXFUNCT
          IF(LUSEBUFF(I))THEN
            DO J=1,NPARAM(I)
              IF(FIX(J,I).EQ.'n')THEN
                NDIM=NDIM+1
                LPAR(J)=XX(NDIM)
              ELSE
                LPAR(J)=PARAM(J,I)
              END IF
            END DO
CADD: introducir la nueva funcion para el calculo de los residuos
c.............................................................................
            IF(CPROF(I).EQ.'c')THEN
              DO K=NF1,NF2
                YP(K,I)=LPAR(1)/
     >           (LPAR(2)*LPAR(2)+(X(K)-LPAR(3))*(X(K)-LPAR(3)))
              END DO
c.............................................................................
            ELSEIF(CPROF(I).EQ.'g')THEN
              DO K=NF1,NF2
                YP(K,I)=LPAR(1)*EXP(-(X(K)-LPAR(3))*(X(K)-LPAR(3))/
     >           (2.*LPAR(2)*LPAR(2)))
              END DO
c.............................................................................
            ELSEIF(CPROF(I).EQ.'d')THEN
              DO K=NF1,NF2
                YP(K,I)=LPAR(1)*EXP(-(X(K)-LPAR(3))*(X(K)-LPAR(3))/
     >           (2.*LPAR(2)*LPAR(2)))
                YP(K,I)=YP(K,I)+FCTEDOUBLET*LPAR(1)*LPAR(2)/LPAR(4)*
     >           EXP(-(X(K)-LPAR(5))*
     >           (X(K)-LPAR(5))/(2.*LPAR(4)*LPAR(4)))
              END DO
c.............................................................................
            ELSEIF(CPROF(I).EQ.'p')THEN
              IF(NPARAM(I).EQ.1)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)
                END DO
              ELSEIF(NPARAM(I).EQ.2)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)+
     >                    LPAR(2)*X(K)
                END DO
              ELSEIF(NPARAM(I).EQ.3)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)+
     >                    LPAR(2)*X(K)+
     >                    LPAR(3)*X(K)*X(K)
                END DO
              ELSEIF(NPARAM(I).EQ.4)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)+
     >                    LPAR(2)*X(K)+
     >                    LPAR(3)*X(K)*X(K)+
     >                    LPAR(4)*X(K)*X(K)*X(K)
                END DO
              ELSEIF(NPARAM(I).EQ.5)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)+
     >                    LPAR(2)*X(K)+
     >                    LPAR(3)*X(K)*X(K)+
     >                    LPAR(4)*X(K)*X(K)*X(K)+
     >                    LPAR(5)*X(K)*X(K)*X(K)*X(K)
                END DO
              ELSEIF(NPARAM(I).EQ.6)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)+
     >                    LPAR(2)*X(K)+
     >                    LPAR(3)*X(K)*X(K)+
     >                    LPAR(4)*X(K)*X(K)*X(K)+
     >                    LPAR(5)*X(K)*X(K)*X(K)*X(K)+
     >                    LPAR(6)*X(K)*X(K)*X(K)*X(K)*X(K)
                END DO
              ELSEIF(NPARAM(I).EQ.7)THEN
                DO K=NF1,NF2
                  YP(K,I)=LPAR(1)+
     >                    LPAR(2)*X(K)+
     >                    LPAR(3)*X(K)*X(K)+
     >                    LPAR(4)*X(K)*X(K)*X(K)+
     >                    LPAR(5)*X(K)*X(K)*X(K)*X(K)+
     >                    LPAR(6)*X(K)*X(K)*X(K)*X(K)*X(K)+
     >                    LPAR(7)*X(K)*X(K)*X(K)*X(K)*X(K)*X(K)
                END DO
              ELSE
                WRITE(*,101)'FATAL ERROR: in subroutine FUNKMFIT.'
                WRITE(*,101)'Polynomial degree too large.'
                STOP
              END IF
c.............................................................................
            END IF
          END IF
        END DO
C
        DFUNK=0.D0
        DO K=NF1,NF2
          FF=0.0
          DO I=1,NMAXFUNCT
            IF(LUSEBUFF(I)) FF=FF+YP(K,I)
          END DO
          DFUNK=DFUNK+DBLE((Y(K)-FF)*(Y(K)-FF))
        END DO
        FUNKMFIT=REAL(DFUNK)
C
101     FORMAT(A)
        END
