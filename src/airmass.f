C------------------------------------------------------------------------------
C Version 28-November-1996                                       file:airmass.f
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
C Program: airmass
C Classification: extinction correction
C Description: Calculate the airmass for fixed observing conditions.
C
Comment
C
        PROGRAM AIRMASSX !note that the program name can not be AIRMASS
                         !since this is already the name of a global variable
                         !in the file redlib.inc
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER READILIM
        REAL READF
C Parameters
        REAL PI
        PARAMETER(PI=3.141593)
C Observatory
        REAL LAT0,LONG0
        REAL LAT,LONG
        REAL ALTOBS
        CHARACTER*15 OBS
C Date
        REAL ANO,MES,DIA,HORA
        REAL TS0
        REAL*8 FJ
        REAL TS0L,TSL
C Exposure time
        REAL HORA0,HORA1,HORA2,DHORA
C Object coordinates, and epoch
        REAL AROBJ0,DECOBJ0
        REAL AROBJ,DECOBJ,TII
        REAL X(3),X0(3)
        REAL ARF,DECF
        REAL ANGH
        REAL ACI,ALT
C Precession elements
        REAL M(3,3)
        REAL GIO,ZETA,TETA
        REAL TFF,TI,TF
        REAL MONTH
C air mass
ccc     REAL AIRMASS
C
        INTEGER K,KK
        INTEGER NINTERV
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL ALT1,ALT2
        REAL HOR0(0:2000),AIRMASS0(0:2000)
        REAL AIRMASS1,AIRMASS2
        REAL AIRMMAX,AIRMMIN,DAIRM
        CHARACTER*10 CHDATE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKFJ1/ANO,MES,DIA,HORA
        COMMON/BLKFJ2/FJ
        COMMON/BLKFJ3/TS0
        COMMON/BLKOBS/LAT,LONG
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        THISPROGRAM='airmass'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        CALL OBSERVAT(OBS,LAT0,LONG0,ALTOBS)
        CALL DECIMAL(LAT0,LAT)
        CALL DECIMAL(LONG0,LONG)
        LONG=LONG/15.
        WRITE(*,101)OBS
        WRITE(*,200)'Longitude: ',LONG
        WRITE(*,200)'Latitude : ',LAT
ccc        WRITE(*,104)
        CALL SHOWHLP('date')
        WRITE(*,101)'DATE'
        WRITE(*,100)'Year           '
        ANO=READF('@')
        WRITE(*,100)'Month (number) '
        MES=READF('@')
        WRITE(*,100)'Day            '
        DIA=READF('@')
        WRITE(CHDATE,'(I2.2,A1,I2.2,A1,I4.4)')INT(MES),'-',INT(DIA),
     +                                        '-',INT(ANO)
        HORA=0.
        CALL FJ_TS0
        WRITE(*,205)'Julian Date: ',FJ
        WRITE(*,200)'TS 0h UT (Greenwich): ',TS0
        TS0L=TS0+LONG
        WRITE(*,200)'TS 0h UT (Local)    : ',TS0L
ccc        WRITE(*,104)
        CALL SHOWHLP('object position')
        WRITE(*,101)'OBJECT POSITION'
        WRITE(*,100)'a.r. (HH.MMSS) '
        AROBJ0=READF('@')
        CALL DECIMAL(AROBJ0,AROBJ)
        WRITE(*,100)'dec. (DD.MMSS) '
        DECOBJ0=READF('@')
        CALL DECIMAL(DECOBJ0,DECOBJ)
        WRITE(*,100)'Equinox '
        TII=READF('@')
C Precession
        IF(INT(MES).LE.2)THEN
          MONTH=AINT((MES-1.)*63./2.)
        ELSE
          MONTH=AINT((MES+1.)*30.6)-63.
        END IF
        MONTH=MONTH+DIA
        TFF=ANO+MONTH/365.25
        TI=(TII-2000.)/100.
        TF=(TFF-2000.-100.*TI)/100.
C
        GIO=((2306.2181+1.39656*TI-.000139*TI*TI)*TF+(.30188-.000344*TI)
     +      *TF*TF+.017998*TF*TF*TF)/3600.
        ZETA=GIO+((.7928+.00041*TI)*TF*TF+.000205*TF*TF*TF)/3600.
        TETA=((2004.3109-.8533*TI-.000217*TI*TI)*TF-(.42665+.000217*TI)
     +       *TF*TF-.041833*TF*TF*TF)/3600.
C
        GIO=GIO*PI/180.
        ZETA=ZETA*PI/180.
        TETA=TETA*PI/180.
C
        M(1,1)=-SIN(GIO)*SIN(ZETA)+COS(GIO)*COS(TETA)*COS(ZETA)
        M(1,2)=-COS(GIO)*SIN(ZETA)-SIN(GIO)*COS(ZETA)*COS(TETA)
        M(1,3)=-SIN(TETA)*COS(ZETA)
        M(2,1)=SIN(GIO)*COS(ZETA)+COS(GIO)*COS(TETA)*SIN(ZETA)
        M(2,2)=COS(GIO)*COS(ZETA)-SIN(GIO)*COS(TETA)*SIN(ZETA)
        M(2,3)=-SIN(TETA)*SIN(ZETA)
        M(3,1)=COS(GIO)*SIN(TETA)
        M(3,2)=-SIN(GIO)*SIN(TETA)
        M(3,3)=COS(TETA)
C
        X0(1)=COS(DECOBJ*PI/180.)*COS(AROBJ*15.*PI/180.)
        X0(2)=COS(DECOBJ*PI/180.)*SIN(AROBJ*15.*PI/180.)
        X0(3)=SIN(DECOBJ*PI/180.)
C
        DO 405,K=1,3
          X(K)=0.
          DO 400,KK=1,3
            X(K)=X(K)+X0(KK)*M(K,KK)
  400     CONTINUE
  405   CONTINUE
        ARF=ATAN2(X(2),X(1))
        DECF=ASIN(X(3))
        ARF=ARF*180./PI/15.
        IF(ARF.LT.0.)ARF=ARF+24.
        DECF=DECF*180./PI
        WRITE(*,*)
        WRITE(*,105)'EPOCH: ',TII
        WRITE(*,200)'a.r.: ',AROBJ
        WRITE(*,200)'dec.: ',DECOBJ
        WRITE(*,*)
        WRITE(*,101)'EPOCH: Observation time'
        WRITE(*,200)'a.r.: ',ARF
        WRITE(*,200)'dec.: ',DECF
C
ccc     WRITE(*,104)
        CALL SHOWHLP('exposure time')
        WRITE(*,101)'EXPOSURE TIME'
        WRITE(*,101)'NOTE: if the final time corresponds to another '//
     +              ' day add 24 hours'
        WRITE(*,100)'UT initial (HH.MMSS) '
        HORA0=READF('@')
        CALL DECIMAL(HORA0,HORA1)
        WRITE(*,100)'UT final   (HH.MMSS) '
        HORA0=READF('@')
        CALL DECIMAL(HORA0,HORA2)
        DHORA=HORA2-HORA1
C
        WRITE(*,100)'No. intervals (max. 2000) '
        NINTERV=READILIM('500',1,2000)
        AIRMASS=0.
        DO K=0,NINTERV
          HOR0(K)=HORA1+DHORA*REAL(K)/REAL(NINTERV)
          TSL=TS0L+HOR0(K)/0.997270
          IF(TSL.GE.24.)TSL=TSL-24.
          ANGH=TSL-ARF
          IF(ANGH.LT.0.)ANGH=ANGH+24.
          CALL CAMBCOOR(ANGH,DECF,ACI,ALT)
          IF(ALT.LE.0.)THEN
            WRITE(*,101)'Object below horizon'
            STOP
          ELSE
            AIRMASS0(K)=1./COS((90.-ALT)*PI/180.)
            AIRMASS=AIRMASS+AIRMASS0(K)
C           WRITE(*,200)'Alture: ',ALT
C           WRITE(*,220)'Time: ',HORA0,',  Air mass: ',AIRMASS
          END IF
          IF(K.EQ.0)THEN
            AIRMASS1=AIRMASS0(K)
            ALT1=ALT
          END IF
          IF(K.EQ.NINTERV)THEN
            AIRMASS2=AIRMASS0(K)
            ALT2=ALT
          END IF
        END DO
        AIRMASS=AIRMASS/REAL(NINTERV+1)
        WRITE(*,240)'Initial values--->  Time: ',HORA1,',     Air mass',
     +              AIRMASS1,'    Alture: ',ALT1
        WRITE(*,240)'Last    values--->  Time: ',HORA2,',     Air mass',
     +              AIRMASS2,'    Alture: ',ALT2
        WRITE(*,230)'Mean air mass: ',AIRMASS
C Plotting
        AIRMMAX=AIRMASS0(0)
        AIRMMIN=AIRMASS0(0)
        DO K=1,NINTERV
          IF(AIRMASS0(K).GT.AIRMMAX)AIRMMAX=AIRMASS0(K)
          IF(AIRMASS0(K).LT.AIRMMIN)AIRMMIN=AIRMASS0(K)
        END DO
        DAIRM=AIRMMAX-AIRMMIN
        AIRMMAX=AIRMMAX+0.05*DAIRM
        AIRMMIN=AIRMMIN-0.05*DAIRM
        HORA2=HORA2+DHORA*0.05
        HORA1=HORA1-DHORA*0.05
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCF(2)
          CALL PGENV(HORA1,HORA2,AIRMMIN,AIRMMAX,0,0)
          CALL PGIDEN_RED
          DO K=0,NINTERV
            CALL PGPOINT(1,HOR0(K),AIRMASS0(K),1)
          END DO
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGMOVE(HORA1,AIRMASS)
          CALL PGDRAW(HORA2,AIRMASS)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGLABEL('Universal Time','Air mass',
     +     OBS//'('//CHDATE//')')
        END DO
        CALL MY_PGEND
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
ccc104     FORMAT(79('='))
105     FORMAT(A,F10.1)
200     FORMAT(A,F10.5)
205     FORMAT(A,F12.2)
ccc220     FORMAT(2(A,F9.5))
230     FORMAT(A,F9.5)
240     FORMAT(3(A,F9.5))
        END

C **********************************************************************
C                                                    SUBROTUINE OBSERVAT
C                                                    *******************
      SUBROUTINE OBSERVAT(OBS,LAT,LONG,ALTOBS)
C
C Longitud y latitud (grados), y altura (metros) de los observatorios
C con el formato DD.MMSS
C Por convenio longitudes son positivas hacia el Este y negativas hacia
C el Oeste.
C
      IMPLICIT NONE
      INCLUDE 'futils.inc'
      INTEGER READI
      REAL READF
C---> argumentos ficticios: SALIDA
      REAL LAT,LONG
      REAL ALTOBS
      CHARACTER*15 OBS
C---> parametros locales
      INTEGER NOBS
      PARAMETER(NOBS=11)
C---> variables locales
      INTEGER I,NOPC
      CHARACTER*15 OBSER(NOBS)
C
      DATA OBSER/'CALAR ALTO','EL TEIDE','MADRID','LA PALMA',
     + 'SAN PEDRO','LICK','LA SILLA','MAUNA KEA','McDonald',
     + 'KPNO','LLANO DEL HATO'/
C
ccc      WRITE(*,104)
      WRITE(*,101)'ENTER OBSERVATORY'
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,160)0,'Observatory data through keyboard (Longitude,'//
     +              'Latitude,Height)'
      DO I=1,NOBS
        IF(I.LE.9)THEN
          WRITE(*,160)I,OBSER(I)
        ELSE
          WRITE(*,161)I,OBSER(I)
        END IF
      END DO
    5 CONTINUE
      CALL SHOWHLP('observatory')
      WRITE(*,*)
      WRITE(*,100) 'Observatory number '
      NOPC=READI('1')
      IF((NOPC.LT.0).OR.(NOPC.GT.NOBS))THEN
        WRITE(*,101)'ERROR: invalid option. Try again.'
        GOTO 5
      END IF
C
      IF(NOPC.EQ.0)THEN
        WRITE(*,*)
        WRITE(*,100)'Observatory name (max. 15 characters)'
        OBS(1:15)=READC('@','@')
        WRITE(*,100)'Latitude (+DD.MMSS)'
        LAT=READF('@')
        WRITE(*,100)'Longitude (+DD.MMSS, +E, -W)'
        LONG=READF('@')
        WRITE(*,100)'Height (metres)'
        ALTOBS=READF('@')
      ELSE IF(NOPC.EQ.1)THEN
C Calar Alto
        LAT=37.1305
        LONG=-2.3240
        ALTOBS=2165.
      ELSE IF(NOPC.EQ.2)THEN
C Izaña
        LAT=28.1732
        LONG=-16.2945
        ALTOBS=2400.
      ELSE IF(NOPC.EQ.3)THEN
C Madrid
        LAT=40.2430
        LONG=-3.4115
        ALTOBS=656.
      ELSE IF(NOPC.EQ.4)THEN
C La Palma
        LAT=28.4540
        LONG=-17.5247
        ALTOBS=2334.
      ELSE IF(NOPC.EQ.5)THEN
C San Pedro (datos aproximados extraidos de un mapa mundi)
        LAT=32.0
        LONG=-115
        ALTOBS=2000.
      ELSE IF(NOPC.EQ.6)THEN
C Lick (datos facilitados por Jesus Gallego)
        LAT=37.2036
        LONG=-123.3812
        ALTOBS=1283.
      ELSE IF(NOPC.EQ.7)THEN
C La Silla
        LAT=-29.1500
        LONG=-70.4400
        ALTOBS=2400.0
      ELSEIF(NOPC.EQ.8)THEN
C UKIRT
        LAT=19.4932
        LONG=-155.2824
        ALTOBS=4194.
      ELSEIF(NOPC.EQ.9)THEN
C McDonald
        LAT=30.4018
        LONG=-104.0118
        ALTOBS=2075.
      ELSEIF(NOPC.EQ.10)THEN
C KPNO
        LAT=31.5748
        LONG=-111.3600
        ALTOBS=2120.
      ELSEIF(NOPC.EQ.11)THEN
C Llano del Hato (Venezuela)
        LAT=08.4724
        LONG=-70.8667
        ALTOBS=3610.
      END IF
      IF(NOPC.NE.0)OBS=OBSER(NOPC)
  100 FORMAT(A,$)
  101 FORMAT(A)
ccc  104 FORMAT(79('='))
  160 FORMAT(1X,'(',I1,')',1X,A)
  161 FORMAT('(',I2,')',1X,A)
      END
C **********************************************************************
C                                                      SUBROUTINE FJ_TS0
C                                                      *****************
      SUBROUTINE FJ_TS0
C
C Calcula la fecha juliana y el Tiempo Sidereo en Greenwich a 0h UT a
C partir de los datos de una fecha concreta. En la fecha juliana se
C tiene en cuenta la fraccion del dia correspondiente a la hora
C traspasada a traves de la variable global HORA.
C Los calculos se realizan en doble precision para no perder informacion
C por redondeo.
C
      IMPLICIT NONE
C---> variables globales: ENTRADA
      REAL ANO,MES,DIA,HORA
C---> variables globales: SALIDA
      REAL TS0
      REAL*8 FJ
C---> variables locales
      INTEGER*4 FECHA
      REAL*8 M,A,B,Y,DT,DTS0
      REAL*8 DANO,DMES,DDIA
C
      COMMON/BLKFJ1/ANO,MES,DIA,HORA
      COMMON/BLKFJ2/FJ
      COMMON/BLKFJ3/TS0
C
      DANO=DBLE(ANO)
      DMES=DBLE(MES)
      DDIA=DBLE(DIA)
      FECHA=INT(ANO)*10000+INT(MES)*100+INT(DIA)
      IF(FECHA.GE.15821015)THEN
        A=DINT(DANO/1.D2)
        B=2.D0-A+DINT(A/4.D0)
      ELSE
        B=0.D0
      END IF
      IF(INT(ANO).GE.0)THEN
        A=0.D0
      ELSE
        A=-.75D0
      END IF
      IF(INT(MES).GE.2)THEN
        Y=DANO
        M=DMES
      ELSE
        Y=DANO-1.D0
        M=DMES+1.2D1
      END IF
      FJ=DINT(365.25D0*Y+A)+DINT(30.6001D0*(M+1.D0))+DDIA+B+
     +   1720994.5D0
      DT=(FJ-2451545.D0)/36525.D0
      DTS0=6.D0*3.6D3+41.D0*6.D1+50.54851D0+8640184.812866D0*DT
     +     +.093104D0*DT*DT-6.2D-6*DT*DT*DT
      DTS0=DMOD(DTS0,8.64D4)
      DTS0=DTS0/8.64D4
      DTS0=DTS0*2.4D1
      TS0=REAL(DTS0)
      IF(TS0.LT.0.)TS0=TS0+24.
      IF(TS0.GT.24.)TS0=TS0-24.
      FJ=FJ+DBLE(HORA/24.)
      END
C **********************************************************************
C                                                    SUBROUTINE CAMBCOOR
C                                                    *******************
      SUBROUTINE CAMBCOOR(HOR,DEC,ACI,ALT)
C
C Transformaciones de coordenadas:
C (angulo horario,declinacion) ----> (acimut,altura)
C
      IMPLICIT NONE
C---> parametros
      REAL PI
      PARAMETER(PI=3.141593)
C---> argumentos ficticios: ENTRADA
      REAL HOR,DEC
C---> argumentos ficticios: SALIDA
      REAL ACI,ALT
C---> variables globales
      REAL LAT,LONG
C---> variables locales
      REAL HORR,DECR
      REAL LATR
      REAL SENOAC,COSEAC
C---> common blocks
      COMMON/BLKOBS/LAT,LONG
C
      HORR=HOR*PI/180.*15.
      DECR=DEC*PI/180.
      LATR=LAT*PI/180.
      SENOAC=COS(DECR)*SIN(HORR)
      COSEAC=-SIN(DECR)*COS(LATR)+COS(DECR)*SIN(LATR)*COS(HORR)
      ACI=ATAN2(SENOAC,COSEAC)
      ACI=ACI*180./PI
      IF(ACI.LT.0.)ACI=ACI+360.
      ALT=ASIN(SIN(LATR)*SIN(DECR)+COS(LATR)*COS(DECR)*COS(HORR))
      ALT=ALT*180./PI
      END
C **********************************************************************
C                                                     SUBROUTINE DECIMAL
C                                                     ******************
      SUBROUTINE DECIMAL(A,B)
C
C Paso de DD.MMSS a DD.dddd
C Utilizamos variables INTEGER*4 para evitar los errores de redondeo.
C Es posible hacer una llamada a la subrutina en la que los argumentos
C verdaderos correspondientes a las argumentos ficticios A y B
C coincidan.
C
      IMPLICIT NONE
C---> argumentos ficticios: ENTRADA
      REAL A
C---> argumentos ficticios: SALIDA
      REAL B
C---> variables locales
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
