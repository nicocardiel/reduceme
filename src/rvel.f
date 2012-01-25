C------------------------------------------------------------------------------
C Version 30-October-1997                                          file: rvel.f
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
C Program: rvel
C Classification: miscellany
C Description: Determines the Earth`s velocity at a fixed time to correct 
C radial velocities.
C
Comment
C
C Calcula la componente de velocidad de la Tierra en una determina direccion
C en la esfera celeste, para los instantes de tiempo solicitados. Para ello 
C hemos utilizamos las formulas descritas en "Astronomical Formulae for 
C Calculators", Jean Meeus.
C
        PROGRAM RVEL
C
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'redlib.inc'
        INTEGER TRUELEN
        INTEGER READILIM
        REAL READF
C
        INTEGER NMAXOBJ,NMAXCOL
        PARAMETER(NMAXOBJ=1000,NMAXCOL=3)
        DOUBLE PRECISION PI
        PARAMETER(PI=3.14159265358979323846D0)
C
        DOUBLE PRECISION F
        DOUBLE PRECISION FJ,JD,JD0,JDEFF
        DOUBLE PRECISION T,L,M,ex
        DOUBLE PRECISION E,KEPLER
        DOUBLE PRECISION NU,LONG,R
        DOUBLE PRECISION CA,CB,CC,CD,CE,CH
        DOUBLE PRECISION EPSILON
        DOUBLE PRECISION X,Y,Z
        DOUBLE PRECISION XX,YY,ZZ
        DOUBLE PRECISION TAU,TAU0
        DOUBLE PRECISION O1,O2,O3
        DOUBLE PRECISION A11,A12,A13,A21,A22,A23,A31,A32,A33
        CHARACTER*11 CFECHA(NMAXCOL),CFECHA_
C
        INTEGER I,II,K,J
        INTEGER NDAYS
        INTEGER NTERMS,NOBJ,NCOL
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL TP(1001),XP(1001),YP(1001),ZP(1001)
        REAL AVX(20),AVY(20),AVZ(20),CHISQR
        REAL XF(1001),YF(1001)
        REAL VELX,VELY,VELZ
        REAL YMIN,YMAX,DY
        REAL RA(NMAXOBJ),DEC(NMAXOBJ),XS,YS,ZS
        REAL RADVEL(NMAXOBJ,NMAXCOL)
        REAL RVOBJ(NMAXOBJ)
        REAL FMEAN(NMAXOBJ)
        CHARACTER*1 CCDIS,CDIS,CPLOT,CFILE,COUT,CREPEAT
        CHARACTER*15 NAME(NMAXOBJ)
        CHARACTER*75 COORFILE,OUTFILE
        LOGICAL LCOLOR(MAX_ID_RED)
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C 
        THISPROGRAM='rvel'
        CALL WELCOME('30-October-1997')
C------------------------------------------------------------------------------
        CREPEAT='n'
        NCOL=0
C
        WRITE(*,100)'Display intermediate data........(y/n) '
        CCDIS(1:1)=READC('n','yn')
        WRITE(*,100)'Plots............................(y/n) '
        CPLOT(1:1)=READC('n','yn')
        IF(CPLOT.EQ.'y')THEN
          CALL PIDEGTER(NTERM,IDN,LCOLOR)
        END IF
C
        F=PI/180.D0
C------------------------------------------------------------------------------
C Entrada de coordenadas de los objectos celestes
        WRITE(*,100)'Enter coordinates from file......(y/n) '
        CFILE(1:1)=READC('y','yn')
        WRITE(*,*)
        IF(CFILE.EQ.'y')THEN
          WRITE(*,101)'* NOTE: File format must be the following:'//
     +     '  FORMAT(F7.4,1X,F8.4,1X,F6.1,1X,A15)'
          WRITE(*,*)
          WRITE(*,101)'      R.A.(J2000.0)DEC  Rvel    Object'
          WRITE(*,101)'      ---------------------------------------'
          WRITE(*,101)'      HH.MMSS +DD.MMSS +VVV.d NameDescription'
          WRITE(*,101)'      123456789012345678901234567890123456789'
          WRITE(*,101)'      000000000111111111122222222223333333333'
          WRITE(*,*)
          WRITE(*,100)'Coordinates file name'
          COORFILE=INFILEX(20,'@',0,0,.0,.0,3,.FALSE.)
          NOBJ=0
20        CONTINUE
          READ(20,'(F7.4,1X,F8.4,1X,F6.1,1X,A15)',END=30)
     +     RA(NOBJ+1),DEC(NOBJ+1),RVOBJ(NOBJ+1),NAME(NOBJ+1)
          NOBJ=NOBJ+1
          IF(NOBJ.EQ.NMAXOBJ)THEN
            WRITE(*,101)'WARNING: NOBJ.EQ.NMAXOBJ'
            GOTO 30
          END IF
          GOTO 20
30        CONTINUE
          CLOSE(20)
        ELSE
          NOBJ=0
40        CONTINUE
          WRITE(*,100)'Object name (max. 15 char.) <ENTER>=EXIT ? '
          READ(*,'(A)')NAME(NOBJ+1)
          IF(TRUELEN(NAME(NOBJ+1)).EQ.0) GOTO 50
          NOBJ=NOBJ+1
          WRITE(*,100)'R.A. (HH.MMSS)  J2000.0 '
          RA(NOBJ)=READF('@')
          WRITE(*,100)'DEC. (DD.MMSS)  J2000.0 '
          DEC(NOBJ)=READF('@')
          WRITE(*,100)'Radial velocity (km/sec)'
          RVOBJ(NOBJ)=READF('@')
          GOTO 40
50        CONTINUE
        END IF
        WRITE(*,*)
        WRITE(*,110)'No. of objects read: ',NOBJ
        WRITE(*,*)
        DO J=1,NOBJ
          CALL DECIMAL(RA(J),RA(J))
          CALL DECIMAL(DEC(J),DEC(J))
          RA(J)=RA(J)*15.*REAL(F)
          DEC(J)=DEC(J)*REAL(F)
        END DO
C------------------------------------------------------------------------------
        WRITE(*,100)'No. of days/10 around date (max. 500) '
        NDAYS=READILIM('30',1,500)
5       WRITE(*,100)'Polynomial degree '
        NTERMS=READILIM('2',2,19)
        IF((NTERMS.LT.2).OR.(NTERMS.GT.19))THEN
          WRITE(*,101)'ERROR: number out of range. Try again.'
          GOTO 5
        END IF
        NTERMS=NTERMS+1
        WRITE(*,*)
C------------------------------------------------------------------------------
10      NCOL=NCOL+1
        WRITE(*,100)'Date.........................(YYYY.MMDDHH) '
        CFECHA_(1:11)=READC('@','1234567890.+-')
        CFECHA(NCOL)=CFECHA_
        JD0=FJ(CFECHA(NCOL))
        JD=FJ('2000.010100')
        IF(CCDIS.EQ.'y')THEN
          WRITE(*,'(A,F17.9)')'FJ0= ',JD0
          WRITE(*,'(A,F17.9)')'FJ = ',JD
        END IF
C------------------------------------------------------------------------------
        DO I=-NDAYS,NDAYS     !bucle para calculo de vector de posicion del Sol
          IF((CCDIS.EQ.'y').AND.(I.EQ.0))THEN
            CDIS='y'
          ELSE
            CDIS='n'
          END IF
          JDEFF=JD0+DBLE(I)/10.D0
          T=(JDEFF-2415020.0D0)/36525.D0
          L=279.69668D0+36000.76892D0*T+.0003025D0*T*T
          M=358.47583D0+35999.04975D0*T-.000150D0*T*T-.0000033D0*T*T*T
          ex=.01675104D0-.0000418D0*T-.000000126D0*T*T
          IF(L.GE.360.D0) L=L-360.D0*DINT(L/360.D0)
          IF(M.GE.360.D0) M=M-360.D0*DINT(M/360.D0)
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.15)')'T  = ',T
            WRITE(*,'(A,F17.13)')'L  = ',L
            WRITE(*,'(A,F17.13)')'M  = ',M
            WRITE(*,'(A,F17.15)')'e  = ',ex
          END IF
          E=KEPLER(M,ex)
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.13)')'E  = ',E
          END IF
          NU=2.D0*DATAN(DSQRT((1.D0+ex)/(1.D0-ex))*DTAN(E*F/2.D0))
          NU=NU/F
          IF(NU.LT.0.D0) NU=NU+360.D0
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.13)')'nu = ',NU
          END IF
          LONG=L+NU-M
          R=1.0000002D0*(1.D0-ex*DCOS(E*F))
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.13)')'lon= ',LONG
            WRITE(*,'(A,F17.15)')'R  = ',R
          END IF
          CA=153.23D0+22518.7541D0*T
          CB=216.57D0+45037.5082D0*T
          CC=312.69D0+32964.3577D0*T
          CD=350.74D0+445267.1142D0*T-0.00144D0*T*T
          CE=231.19D0+20.20D0*T
          CH=353.40D0+65928.7155D0*T
          CA=CA*F
          CB=CB*F
          CC=CC*F
          CD=CD*F
          CE=CE*F
          CH=CH*F
          LONG=LONG+0.00134D0*DCOS(CA)+
     +              0.00154D0*DCOS(CB)+
     +              0.00200D0*DCOS(CC)+
     +              0.00179D0*DSIN(CD)+
     +              0.00178D0*DSIN(CE)
          R=R+0.00000543D0*DSIN(CA)+
     +        0.00001575D0*DSIN(CB)+
     +        0.00001627D0*DSIN(CC)+
     +        0.00003076D0*DCOS(CD)+
     +        0.00000927D0*DSIN(CH)
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.13)')'lon= ',LONG
            WRITE(*,'(A,F17.15)')'R  = ',R
          END IF
          EPSILON=23.452294D0-0.0130125D0*T-0.00000164D0*T*T+
     +            0.000000503D0*T*T*T
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.13)')'eps= ',EPSILON
          END IF
          X=R*DCOS(LONG*F)
          Y=R*DSIN(LONG*F)*DCOS(EPSILON*F)
          Z=R*DSIN(LONG*F)*DSIN(EPSILON*F)
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.15)')'X  = ',X
            WRITE(*,'(A,F17.15)')'Y  = ',Y
            WRITE(*,'(A,F17.15)')'Z  = ',Z
          END IF
          TAU0=(JDEFF-2451544.5333981D0)/36525.D0    !usamos formulas recientes
          TAU=(JD-JDEFF)/36525.D0
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.14)')'t0 = ',TAU0
            WRITE(*,'(A,F17.14)')'t  = ',TAU
          END IF
          CALL PRECE(TAU0,TAU,O1,O2,O3)
          O1=O1/3600.D0
          O2=O2/3600.D0
          O3=O3/3600.D0
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.15)')'O1 = ',O1
            WRITE(*,'(A,F17.15)')'O2 = ',O2
            WRITE(*,'(A,F17.15)')'O3 = ',O3
          END IF
          O1=O1*F
          O2=O2*F
          O3=O3*F
          A11=DCOS(O1)*DCOS(O2)*DCOS(O3)-DSIN(O1)*DSIN(O2)
          A21=DSIN(O1)*DCOS(O2)+DCOS(O1)*DSIN(O2)*DCOS(O3)
          A31=DCOS(O1)*DSIN(O3)
          A12=-DCOS(O1)*DSIN(O2)-DSIN(O1)*DCOS(O2)*DCOS(O3)
          A22=DCOS(O1)*DCOS(O2)-DSIN(O1)*DSIN(O2)*DCOS(O3)
          A32=-DSIN(O1)*DSIN(O3)
          A13=-DCOS(O2)*DSIN(O3)
          A23=-DSIN(O2)*DSIN(O3)
          A33=DCOS(O3)
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.15)')'A11= ',A11
            WRITE(*,'(A,F17.15)')'A12= ',A12
            WRITE(*,'(A,F17.15)')'A13= ',A13
            WRITE(*,'(A,F17.15)')'A21= ',A21
            WRITE(*,'(A,F17.15)')'A22= ',A22
            WRITE(*,'(A,F17.15)')'A23= ',A23
            WRITE(*,'(A,F17.15)')'A31= ',A31
            WRITE(*,'(A,F17.15)')'A32= ',A32
            WRITE(*,'(A,F17.15)')'A33= ',A33
          END IF
          XX=A11*X+A12*Y+A13*Z
          YY=A21*X+A22*Y+A23*Z
          ZZ=A31*X+A32*Y+A33*Z
C cambio de origen de la Tierra al Sol
          XX=-XX
          YY=-YY
          ZZ=-ZZ
          IF(CDIS.EQ.'y')THEN
            WRITE(*,'(A,F17.9)')'FJ = ',JDEFF
            WRITE(*,'(A,F17.15)')'X  = ',XX
            WRITE(*,'(A,F17.15)')'Y  = ',YY
            WRITE(*,'(A,F17.15)')'Z  = ',ZZ
          END IF
          XP(I+NDAYS+1)=REAL(XX)
          YP(I+NDAYS+1)=REAL(YY)
          ZP(I+NDAYS+1)=REAL(ZZ)
          TP(I+NDAYS+1)=REAL(I)
        END DO                                                   !fin del bucle
C------------------------------------------------------------------------------
C Dibujamos posiciones calculadas
        IF(CPLOT.EQ.'y')THEN
          YMIN=XP(1)
          YMAX=YMIN
          DO I=2,2*NDAYS+1
            IF(XP(I).LT.YMIN) YMIN=XP(I)
            IF(XP(I).GT.YMAX) YMAX=XP(I)
            IF(YP(I).LT.YMIN) YMIN=YP(I)
            IF(YP(I).GT.YMAX) YMAX=YP(I)
            IF(ZP(I).LT.YMIN) YMIN=ZP(I)
            IF(ZP(I).GT.YMAX) YMAX=ZP(I)
          END DO
          DY=YMAX-YMIN
          DY=DY/50.
          YMIN=YMIN-DY*3.
          YMAX=YMAX+DY*3.
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSCF(2)
            CALL PGENV(-REAL(NDAYS+1),REAL(NDAYS+1),YMIN,YMAX,0,0)
            CALL PGLABEL('days/10','X,Y,Z','Earth position ['//
     +       CFECHA(NCOL)//']')
            CALL PGPOINT(2*NDAYS+1,TP,XP,17)
            CALL PGSCH(1.5)
            CALL PGPTEXT(TP(1),XP(1)+DY,0.,.5,'X')
            CALL PGSCH(1.0)
            CALL PGPOINT(2*NDAYS+1,TP,YP,17)
            CALL PGSCH(1.5)
            CALL PGPTEXT(TP(1),YP(1)+DY,0.,.5,'Y')
            CALL PGSCH(1.0)
            CALL PGPOINT(2*NDAYS+1,TP,ZP,17)
            CALL PGSCH(1.5)
            CALL PGPTEXT(TP(1),ZP(1)+DY,0.,.5,'Z')
            CALL PGSCH(1.0)
          END DO
        END IF
C------------------------------------------------------------------------------
C Ajustamos polinomios a las posiciones calculadas
        CALL POLFIT(TP,XP,XP,2*NDAYS+1,NTERMS,0,AVX,CHISQR)
        CALL POLFIT(TP,YP,YP,2*NDAYS+1,NTERMS,0,AVY,CHISQR)
        CALL POLFIT(TP,ZP,ZP,2*NDAYS+1,NTERMS,0,AVZ,CHISQR)
C------------------------------------------------------------------------------
C Dibujamos polinomios ajustados
        IF(CPLOT.EQ.'y')THEN
          DO II=-NDAYS,NDAYS
            I=II+NDAYS+1
            XF(I)=REAL(II)
            YF(I)=AVX(NTERMS)
            DO K=NTERMS-1,1,-1
              YF(I)=YF(I)*XF(I)+AVX(K)
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGLINE(2*NDAYS+1,XF,YF)
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
          END DO
          DO II=-NDAYS,NDAYS
            I=II+NDAYS+1
            XF(I)=REAL(II)
            YF(I)=AVY(NTERMS)
            DO K=NTERMS-1,1,-1
              YF(I)=YF(I)*XF(I)+AVY(K)
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGLINE(2*NDAYS+1,XF,YF)
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
          END DO
          DO II=-NDAYS,NDAYS
            I=II+NDAYS+1
            XF(I)=REAL(II)
            YF(I)=AVZ(NTERMS)
            DO K=NTERMS-1,1,-1
              YF(I)=YF(I)*XF(I)+AVZ(K)
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGLINE(2*NDAYS+1,XF,YF)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END IF
C------------------------------------------------------------------------------
C Calculamos las derivadas de los polinomios
ccc60      I=NDAYS+1                               !calculamos en el punto central
        I=NDAYS+1                               !calculamos en el punto central
        XF(I)=0.
        VELX=AVX(NTERMS)*REAL(NTERMS-1)
        VELY=AVY(NTERMS)*REAL(NTERMS-1)
        VELZ=AVZ(NTERMS)*REAL(NTERMS-1)
        DO K=NTERMS-1,2,-1
          VELX=VELX*XF(I)+AVX(K)*REAL(K-1)
          VELY=VELY*XF(I)+AVY(K)*REAL(K-1)
          VELZ=VELZ*XF(I)+AVZ(K)*REAL(K-1)
        END DO
        VELX=VELX*10.
        VELY=VELY*10.
        VELZ=VELZ*10.
        IF(CCDIS.EQ.'y')THEN
          WRITE(*,'(A,F9.7,A)')'Vx = ',VELX,' AU/day'
          WRITE(*,'(A,F9.7,A)')'Vy = ',VELY,' AU/day'
          WRITE(*,'(A,F9.7,A)')'Vz = ',VELZ,' AU/day'
        END IF
C Pasamos de AU/day a km/s  (1 UA=1.49597870E11 metros)
        VELX=VELX*1.49597870E8/86400.
        VELY=VELY*1.49597870E8/86400.
        VELZ=VELZ*1.49597870E8/86400.
        IF(COUT.EQ.'y') WRITE(30,101)'Date: '//CFECHA(NCOL)
        DO J=1,NOBJ
          XS=COS(RA(J))*COS(DEC(J))
          YS=SIN(RA(J))*COS(DEC(J))
          ZS=SIN(DEC(J))
          RADVEL(J,NCOL)=XS*VELX+YS*VELY+ZS*VELZ
          IF(CCDIS.EQ.'y')THEN
            WRITE(*,'(A,F6.2,A)')'Object: '//NAME(J)//'   '//
     +       '>>> Velocity = ',RADVEL(J,NCOL),' (km/sec)'
          END IF
        END DO
C------------------------------------------------------------------------------
        IF(NCOL.LT.NMAXCOL)THEN
          WRITE(*,100)'New date (y/n) '
          CREPEAT(1:1)=READC('y','yn')
          IF(CREPEAT.EQ.'y') GOTO 10
        END IF
C------------------------------------------------------------------------------
C Generamos fichero de salida con toda la informacion
        IF(CCDIS.EQ.'y')THEN
          WRITE(*,100)'Create output file name (y/n) '
          COUT(1:1)=READC('n','yn')
        ELSE
          COUT='y'
        END IF
        IF(COUT.EQ.'y')THEN
          WRITE(*,*)
          WRITE(*,100)'Output file name'
          OUTFILE=OUTFILEX(30,'@',0,0,.0,.0,3,.FALSE.)
          WRITE(30,101)'File: '//OUTFILE(1:TRUELEN(OUTFILE))//
     +     '  [Output from program: rvel]'
          IF(CFILE.EQ.'y')THEN
            WRITE(30,101)'* File name with coordinates: '//
     +       COORFILE(1:TRUELEN(COORFILE))
          ELSE
            WRITE(30,101)'* NOTE: object coordinates have been '//
     +       'entered through keyboard.'
          END IF
          WRITE(30,101)'* NOTE: radial velocities are given in km/sec'
          WRITE(30,*)
          WRITE(30,100)'Object            RadVel'
          DO K=1,NCOL
            WRITE(30,'(4X,A11,$)')CFECHA(K)
          END DO
          WRITE(30,'(4X,A4)')'mean'
          WRITE(30,101)'---------------------------------------'//
     +     '----------------------------------------'
          DO J=1,NOBJ
            WRITE(30,'(A15,1X,F8.1,$)')NAME(J),RVOBJ(J)
            FMEAN(J)=0.
            DO K=1,NCOL 
              WRITE(30,'(1X,F8.1,1X,F5.1,$)')
     +         RVOBJ(J)-RADVEL(J,K),RADVEL(J,K)
              FMEAN(J)=FMEAN(J)+(RVOBJ(J)-RADVEL(J,K))
            END DO
            FMEAN(J)=FMEAN(J)/REAL(NCOL)
            WRITE(30,'(1X,F8.1)') FMEAN(J)
          END DO
          CLOSE(30)
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Create output file for index (y/n) '
        COUT(1:1)=READC('n','yn')
        IF(COUT.EQ.'y')THEN
          WRITE(*,100)'File name'
          OUTFILE=OUTFILEX(32,'@',0,0,.0,.0,3,.FALSE.)
          DO J=1,NOBJ
            WRITE(32,'(F8.1,1X,A15)')FMEAN(J),NAME(J)
          END DO
          CLOSE(32)
        END IF
C------------------------------------------------------------------------------
        IF(CPLOT.EQ.'y') CALL PGEND
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Transforma de DD.MMSS a DD.dddd
        SUBROUTINE DECIMAL(A,B)
        IMPLICIT NONE
        REAL A,B
        INTEGER*4 AA,A1,A2,A3
        CHARACTER*1 SIGNO
C
        SIGNO='+'
        IF(A.LT.0)SIGNO='-'
        AA=ABS(NINT(A*10000))
        A1=AA/10000
        A2=AA-A1*10000
        A2=A2/100
        A3=AA-A1*10000-A2*100
        B=REAL(A1)+REAL(A2)/60.+REAL(A3)/3600.
        IF(SIGNO.EQ.'-')B=-B
        END

C
C******************************************************************************
C
        DOUBLE PRECISION FUNCTION FJ(CFECHA)
C
        IMPLICIT NONE
        CHARACTER*11 CFECHA
C
        INTEGER IANO,IMES,IDIA,IHORA
        INTEGER FECHA
        DOUBLE PRECISION ANO,MES,DIA,HORA
        DOUBLE PRECISION A,B,Y,M
C------------------------------------------------------------------------------
        READ(CFECHA,'(I4,1X,I2,I2,I2)')IANO,IMES,IDIA,IHORA
        ANO=DBLE(IANO)
        MES=DBLE(IMES)
        DIA=DBLE(IDIA)
        HORA=DBLE(IHORA)
        FECHA=IANO*10000+IMES*100+IDIA
        IF(FECHA.GE.15821015)THEN
          A=DINT(ANO/1.D2)
          B=2.D0-A+DINT(A/4.D0)
        ELSE
          B=0.D0
        END IF
        IF(IANO.GE.0)THEN
          A=0.D0
        ELSE
          A=-.75D0
        END IF
        IF(IMES.GE.2)THEN
          Y=ANO
          M=MES
        ELSE
          Y=ANO-1.D0
          M=MES+1.2D1
        END IF
        FJ=DINT(365.25D0*Y+A)+DINT(30.6001D0*(M+1.D0))+DIA+B+
     +     1720994.5D0
        FJ=FJ+HORA/24.D0
C
        END
C
C******************************************************************************
C
        DOUBLE PRECISION FUNCTION KEPLER(M,ex)
        IMPLICIT NONE
        DOUBLE PRECISION PI
        PARAMETER(PI=3.14159265358979323846D0)
C
        INTEGER N
        DOUBLE PRECISION M,ex
        DOUBLE PRECISION E0,E1,F
C------------------------------------------------------------------------------
        F=PI/180.D0
        E0=M
        N=0
10      CONTINUE
        N=N+1
        E1=E0+(M+ex/F*DSIN(E0*F)-E0)/(1.D0-ex*DCOS(E0*F))
        IF(N.GT.5000) GOTO 20
        IF(DABS(E1-E0).GT.1.0D-9)THEN
          E0=E1
          GOTO 10
        END IF
20      CONTINUE
        KEPLER=E1
        END
C
C******************************************************************************
C Calcula los angulos de rotacion por efecto de precesion
C (formulas extraidas de "Computatinal Spherical Astronomy", Laurence G. Taff,
C pag. 35
        SUBROUTINE PRECE(T1,T2,O1,O2,O3)
        IMPLICIT NONE
        DOUBLE PRECISION T1,T2
        DOUBLE PRECISION O1,O2,O3
C------------------------------------------------------------------------------
        O1=(2306.2181D0+1.39656D0*T1-0.000139D0*T1*T1)*T2+
     +     (0.30188D0-0.000344D0*T1)*T2*T2+0.017998D0*T2*T2*T2
        O2=O1+(0.79280D0+0.000410D0*T1)*T2*T2+0.000205D0*T2*T2*T2
        O3=(2004.3109D0-0.85330D0*T1-0.000217D0*T1*T1)*T2
     +     -(0.42665D0+0.000217D0*T1)*T2*T2-0.041833D0*T2*T2*T2
        END
C
C Version 23-Mayo-1995
C******************************************************************************
C                                              INTEGER FUNCTION TRUELEN(CADENA)
C******************************************************************************
C
        INTEGER FUNCTION TRUELEN(CADENA)
C busca la longitud de una cadena ignorando espacio en blanco al final
        IMPLICIT NONE
        CHARACTER*(*) CADENA
        INTEGER I,L
C
        L=LEN(CADENA)
C
        DO I=L,1,-1
          IF(ICHAR(CADENA(I:I)).NE.32)THEN
            TRUELEN=I
            RETURN
          END IF
        END DO
        TRUELEN=0
        END
