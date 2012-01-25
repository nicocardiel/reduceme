C------------------------------------------------------------------------------
C Version 28-November-1996                                     file: chancoor.f
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
C Program: chancoor
C Classification: miscellany
C Description: Transform r.a. and dec. coordinates from an initial equinox 
C to another equinox.
C
Comment
C
        PROGRAM CHANCOOR
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        REAL PI
        PARAMETER (PI=3.141593)
        INTEGER I,K,KK
        REAL TI,TF,TII,TFF
        REAL GIO,ZETA,TETA
        REAL X(3),X0(3)
        REAL M(3,3)
        REAL ARH,ARM,ARS
        REAL DECD,DECM,DECS
        CHARACTER*1 SIGNO
        CHARACTER*1 CCONT
!       CHARACTER*80 CLINEA
        REAL ARI,ARF
        REAL DECI,DECF
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='chancoor'
        CALL WELCOME('6-December-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Initial Equinox'
        TII=READF('@')
        WRITE(*,100)'Final   Equinox'
        TFF=READF('@')
        TI=(TII-2000.)/100.
        TF=(TFF-2000.-100.*TI)/100.
C elementos precesionales en grados
        GIO=((2306.2181+1.39656*TI-.000139*TI*TI)*TF+(.30188-.000344*TI)
     +    *TF*TF+.017998*TF*TF*TF)/3600.
        ZETA=GIO+((.7928+.00041*TI)*TF*TF+.000205*TF*TF*TF)/3600.
        TETA=((2004.3109-.8533*TI-.000217*TI*TI)*TF-(.42665+.000217*TI)
     +     *TF*TF-.041833*TF*TF*TF)/3600.
C los pasamos a radianes
        GIO=GIO*PI/180.
        ZETA=ZETA*PI/180.
        TETA=TETA*PI/180.
C matriz de rotacion
        M(1,1)=-SIN(GIO)*SIN(ZETA)+COS(GIO)*COS(TETA)*COS(ZETA)
        M(1,2)=-COS(GIO)*SIN(ZETA)-SIN(GIO)*COS(ZETA)*COS(TETA)
        M(1,3)=-SIN(TETA)*COS(ZETA)
        M(2,1)=SIN(GIO)*COS(ZETA)+COS(GIO)*COS(TETA)*SIN(ZETA)
        M(2,2)=COS(GIO)*COS(ZETA)-SIN(GIO)*COS(TETA)*SIN(ZETA)
        M(2,3)=-SIN(TETA)*SIN(ZETA)
        M(3,1)=COS(GIO)*SIN(TETA)
        M(3,2)=-SIN(GIO)*SIN(TETA)
        M(3,3)=COS(TETA)
        DO I=1,3
          WRITE(*,*)M(I,1),M(I,2),M(I,3)
        END DO
C entramos coordenadas del objeto
10      WRITE(*,*)
        WRITE(*,101)'*******************************'
        WRITE(*,100)'A.R. (HH,MM,SS)? '
!       CLINEA=READC('@','01234567890,')
!       READ(CLINEA,*)ARH,ARM,ARS
        READ(*,*)ARH,ARM,ARS
        ARI=ARH+ARM/60+ARS/3600
        WRITE(*,100)'DEC. SIGN (+/-)'
        SIGNO(1:1)=READC('+','+-')
        WRITE(*,100)'DEC. (DD,MM,SS)? '
!       CLINEA=READC('@','01234567890,')
!       READ(CLINEA,*)DECD,DECM,DECS
        READ(*,*)DECD,DECM,DECS
        DECI=DECD+DECM/60+DECS/3600
        IF(SIGNO.EQ.'-')DECI=-DECI
C coordenadas rectangulares del objeto
        X0(1)=COS(DECI*PI/180.)*COS(ARI*15.*PI/180.)
        X0(2)=COS(DECI*PI/180.)*SIN(ARI*15.*PI/180.)
        X0(3)=SIN(DECI*PI/180.)
C cambio a coordenadas de la epoca
        DO K=1,3
          X(K)=0.
          DO KK=1,3
            X(K)=X(K)+X0(KK)*M(K,KK)
          END DO
        END DO
        ARF=ATAN2(X(2),X(1))
        DECF=ASIN(X(3))
        ARF=ARF*180./PI/15.
        IF(ARF.LT.0.)ARF=ARF+24.
        DECF=DECF*180./PI
C
        ARH=AINT(ARF)
        ARM=(ARF-ARH)*60
        ARS=(ARM-AINT(ARM))*60
        ARM=AINT(ARM)
C
        IF(DECF.LT.0)THEN
          SIGNO='-'
          DECF=-DECF
        ELSE
          SIGNO='+'
        END IF
        DECD=AINT(DECF)
        DECM=(DECF-DECD)*60
        DECS=(DECM-AINT(DECM))*60
        DECM=AINT(DECM)
C
        WRITE(*,117)'AR (HH MM SS): ',INT(ARH),INT(ARM),ARS
        WRITE(*,111)'DEC(DD MM SS): ',SIGNO,INT(DECD),INT(DECM),DECS
        WRITE(*,100)'Continue (y/n) '
        CCONT(1:1)=READC('y','yn')
        IF((CCONT.EQ.'N').OR.(CCONT.EQ.'n'))THEN
        ELSE
          GOTO 10
        END IF
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
117     FORMAT(A,3X,I2.2,1X,I2.2,1X,F4.1)
111     FORMAT(A,2X,A1,I2.2,1X,I2.2,1X,F4.1)
        END
