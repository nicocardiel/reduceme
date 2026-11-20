!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: loglin.f
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
! Program: loglin
! Classification: miscellany
! Description: Rebins data from a linear wavelength coverage to logarithmic 
! scale
!
!Comment
!
! Convierte un espectro de escala lineal a escala logaritmica
!
        PROGRAM LOGLIN10
!
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
!
        INTEGER I,J
        INTEGER MIN,MAX,TLOG
        REAL S(NCMAX)
        REAL STWV2,DISP2
        CHARACTER*255 INFILE,OUTFILE
!
        COMMON/BLK1/S
        COMMON/BLK2/MIN,MAX
        COMMON/BLK3/NCHAN,STWV,DISP,STWV2,DISP2
!------------------------------------------------------------------------------
        THISPROGRAM='loglin'
        CALL WELCOME('20-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Input  file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
!
        WRITE(*,100,ADVANCE='NO') '1=NATURAL LOG , 2=BASE 10 LOG'
        TLOG=READILIM('@',1,2)

        WRITE(*,100,ADVANCE='NO')'Minimum channel with data '
        MIN=READILIM('@',1,NCHAN)
        WRITE(*,100,ADVANCE='NO')'Maximum channel with data '
        MAX=READILIM('@',MIN,NCHAN)
!
        call bin(0,TLOG)
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV2,DISP2,1,.FALSE.)

        DO I=1,NSCAN
          READ(20) (S(J),J=1,NCHAN)
          CALL BIN(1,TLOG)
          WRITE(30) (S(J),J=1,NCHAN)
        END DO
!
        CLOSE(20)
        CLOSE(30)
!
        STOP
100     FORMAT(A)
        END
!
!******************************************************************************
!
        SUBROUTINE BIN(MODE,TLOG)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER MJ,I,MODE
        INTEGER MIN,MAX,TLOG
        INTEGER LMIN,LMIN2,LMIN3
        REAL S(NCMAX)
        REAL WMIN
        REAL STWV2,DISP2
        DOUBLE PRECISION W(NCMAX+1)
        DOUBLE PRECISION S2(NCMAX)       !espectro de salida (esc. logaritmica)
        DOUBLE PRECISION S3(NCMAX)           !espectro de entrada (esc. lineal)
        DOUBLE PRECISION DLOGW,DW,DV,WL1,WL2,EXC,RES
        DOUBLE PRECISION DL1,DL2,XCNT1,XCNT2,DIFF
!
        COMMON/BLK1/S
        COMMON/BLK2/MIN,MAX
        COMMON/BLK3/NCHAN,STWV,DISP,STWV2,DISP2
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        WMIN=STWV-DISP/2.              !l.d.o. borde izquierdo del primer pixel

        IF(MODE.EQ.1) THEN
          DO MJ=1,NCHAN
            S3(MJ)=DBLE(S(MJ))
          END DO
        END IF
        DW=DBLE(DISP)
        DO MJ=1,NCHAN+1          !modificado en 1 respecto a la rutina original
          W(MJ)=DBLE(WMIN)+DBLE(MJ-1)*DW    !l.d.o. borde izquierdo del pixel J
        END DO
! En la subrutina en el Vax la siguiente linea fallaba si MAX=NCHAN, pero al
! sobredimensionar W en 1, ahora existe W(NCMAX+1) y la instruccion es valida
! tambien en ese caso.
        IF(TLOG.EQ.1) THEN
          DLOGW=(DLOG(W(MAX+1)/W(MIN)))/DBLE(NCHAN)             !wavelength bin
        ELSE
          DLOGW=(DLOG10(W(MAX+1)/W(MIN)))/DBLE(NCHAN)           !wavelength bin
        END IF
        DISP2=REAL(DLOGW)
        DV=(299793.D0)*DLOGW
        IF(TLOG.EQ.2) DV=DV/dble(log10(exp(1.)))         !velocity bin (km/sec)
        IF(MODE.EQ.0) THEN
          WRITE(*,100,ADVANCE='NO')'>>> Wavelength bin.......: '
          WRITE(*,*)DLOGW
          WRITE(*,100,ADVANCE='NO')'>>> Velocity bin (km/sec): '
          WRITE(*,*)DV
        END IF
!
        IF(TLOG.EQ.1) THEN
          DO MJ=1,NCHAN+1
            W(MJ)=DLOG(W(MJ))   !log. neperianos l.d.o. borde izquierdo pixel J
          END DO
        ELSE
          DO MJ=1,NCHAN+1
            W(MJ)=DLOG10(W(MJ))    !log. decimal l.d.o. borde izquierdo pixel J
          END DO
        END IF
        STWV2=W(MIN)+DISP2/2.

        IF(MODE.EQ.0) RETURN
!
        DO MJ=1,NCHAN
          S2(MJ)=0.0D+00
        END DO
!
! BIN TO LOG WAVELENGTH
!
! RES es el numero de cuentas del pixel anterior (lineal) no usadas en la
! generacion del ultimo pixel (logaritmico)
        RES=0.0D+00
        WL1=W(MIN)                         !LOG(l.d.o. inicial del nuevo pixel)
        WL2=WL1+DLOGW                        !LOG(l.d.o. final del nuevo pixel)
        LMIN=MIN
        LMIN2=LMIN+1
        IF(WL2.GE.W(LMIN2)) RES=S3(MIN)
        DO 12 MJ=1,NCHAN
          LMIN3=LMIN2+1
! si el nuevo pixel (logaritmico) contiene mas de un antiguo pixel (lineal)
! se van sumando. Se suma hasta el ultimo pixel (lineal) completo abarcado,
! es decir, siempre se suma por defecto
          DO 13 I=LMIN3,MAX
            LMIN2=I-1
            IF(W(I).GT.WL2) GO TO 14
            S2(MJ)=S3(LMIN2)+S2(MJ)
 13       CONTINUE
 14       CONTINUE
          S2(MJ)=S2(MJ)+RES
          IF(W(LMIN2).GE.WL2) GO TO 15
          DL1=W(LMIN2+1)-W(LMIN2)
          DL2=WL2-W(LMIN2)
          EXC=S3(LMIN2)*(DL2/DL1)
          RES=S3(LMIN2)-EXC
          S2(MJ)=S2(MJ)+EXC
          WL1=WL2
          WL2=WL2+DLOGW
          IF(WL2.LE.W(LMIN2+1)) RES=0.0D+00
          LMIN=LMIN2
          LMIN2=LMIN+1
          GO TO 20
 15       DL1=W(LMIN2)-W(LMIN)
          DL2=WL2-WL1
          S2(MJ)=S2(MJ)+S3(LMIN)*(DL2/DL1)
          WL1=WL2
          WL2=WL2+DLOGW
          IF(WL2.GT.W(LMIN2)) GO TO 19
          RES=0.0D+00
          GO TO 20
 19       DL1=W(LMIN2)-W(LMIN)
          DL2=W(LMIN2)-WL1
          RES=S3(LMIN)*(DL2/DL1)
 20       CONTINUE
 12     CONTINUE
!
! ADD UP COUNTS TO CHECK
!
        XCNT1=0.0D0
        XCNT2=0.0D0
        DO I=MIN,MAX
          XCNT1=XCNT1+S3(I)
        END DO
        DO I=1,NCHAN
          XCNT2=XCNT2+S2(I)
        END DO
        WRITE(*,100,ADVANCE='NO')'Check number of counts: '
        DIFF=DABS(XCNT2-XCNT1)
        WRITE(*,*) XCNT1,XCNT2,DIFF
!
        DO I=1,NCHAN
          S(I)=S2(I)
        END DO
!
100     FORMAT(A)
        END
