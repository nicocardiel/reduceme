C------------------------------------------------------------------------------
C Version 13-January-2009                                       File: mideind.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C INTEGER FUNCTION MIDEIND(NS1,NS2,ITI,WV,FWV,CERR,RCVEL1,NCRES,FINDEX,EINDEX,
C                          EJJGG,ESIMU,SN,FFPLAW,LONLY_SN)
C
C Input: NS1,NS2,ITI,WV,CERR,RCVEL1,NCRES,FFPLAW,LONLY_SN
C Output: FINDEX,EINDEX,EJJGG,ESIMU,SN
C Input/Output: (see COMMON blocks)
C
C Return MIDEIND=0 if atomic/molecular/D4000/B4000/generic index (and error)
C has been properly measured. MIDEIND=1 otherwise.
C
C INTEGER     NS1,NS2 -> scans to be coadded
C INTEGER     ITI -> index type:
C                    ITI= -??: slope index
C                              ITI = -C, where C: number of continuum regions 
C                              Since Cmin=2 and Cmax=99
C                              ==> ITImin=-2 and ITImax=-99
C                    ITI=   1: molecular index
C                    ITI=   2: atomic index
C                    ITI=   3: D4000 (Bruzual 1983)
C                    ITI=   4: B4000 (own defintion)
C                    ITI =  5: like B4000 but computing flux per angstrom
C                    ITI=????: generic index:
C                              ITI= C x 100 + L, where
C                                   C: number of continuum regions
C                                   L: number of absorption regions
C                              Since Cmin=1, Cmax=99, Lmin=1, Lmax=99
C                              ==> ITImin=101, ITImax=9999
C REAL        WV(NWVMAX) -> wavelength limits
C REAL        FWV(NWVMAX/4) -> constant factors to be used when computing
C                           the index as multiplicative coefficients for
C                           the absorption signal.
C CHARACTER*1 CERR -> 'y' if error spectrum is available, 'n' otherwise
C REAL        RCVEL1 -> 1 + v/c
C INTEGER     NCRES -> No. of response curve to be employed (1=averaged)
C REAL        FINDEX -> measured index
C REAL        EINDEX -> measured index error using own formulae
C REAL        EJJGG -> measured index error using JJGGs formulae
C REAL        ESIMU -> measured index error using numerical simulations
C REAL        SN -> averaged (S/N ratio)/Angstrom in the index region
C REAL        FFPLAW -> =0,...,10 indicates the fraction of light (in the
C                    selected photometric band) to be used with a power law
C                    of the form: F(lambda) = k * (lambda)**(alpha-2.0),
C                    where alpha is defined through the variable ALPHAPLAW
C                    (global variable ---COMMON---). If FFPLAW = -1 no
C                    power law is employed.
C LOGICAL     LONLY_SN -> if .TRUE., the subroutine only computes the
C                         mean S/N ratio and returns
C
Comment
C------------------------------------------------------------------------------
        INTEGER FUNCTION MIDEIND(NS1,NS2,ITI,WV,FWV,CERR,RCVEL1,NCRES,
     +   FINDEX,EINDEX,EJJGG,ESIMU,SN,FFPLAW,LONLY_SN)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        INTEGER NFMAX                        !numero maximo de curvas respuesta
        PARAMETER (NFMAX=101)
C
C parametros de la subrutina:
        INTEGER NS1,NS2                               !numeros de scans a sumar
        INTEGER ITI                                             !tipo de indice
        REAL WV(NWVMAX)                                   !l.d.o. de las bandas
        REAL FWV(NWVMAX/4)   !constantes para multiplicar se\~{n}al de asbsorc.
        CHARACTER*1 CERR                                    !calculo de errores 
        REAL RCVEL1                                         !.............(1+z)
        INTEGER NCRES                !numero de la curva respuesta (1=promedio)
        REAL FINDEX,EINDEX                                      !indice y error
        REAL EJJGG           !error calculado con las formulas de JJGG's Thesis
        REAL ESIMU      !error que se obtiene al simular numericamente el error
        REAL SN   !valor medio de la relacion senhal/ruido en las bandas usadas
        INTEGER FFPLAW        !indica la fraccion de luz de la ley de potencias
        LOGICAL LONLY_SN          !indica si se desea solamente la relacion S/N
C
C variables globales (visibles a traves de COMMONs)
        INTEGER NCREST
        INTEGER NSIMUL,NSEED
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL FLUX(NCMAX,NFMAX)
        REAL WLMIN
        REAL ALPHAPLAW,EBVPLAW
        REAL RESCHAN(NCMAX)
        CHARACTER*1 CBANDPLAW
        LOGICAL LPLOT,LPERRI
        LOGICAL LCOLOR(MAX_ID_RED)
C
C parametros
        REAL CTE
        PARAMETER(CTE=1.085736205)                         !2.5*ALOG10(EXP(1.))
C
C variables locales:
        INTEGER I,J,JJ,K
        INTEGER NBAND,NB,NBB
        INTEGER J1(NBDMAX),J2(NBDMAX)
        INTEGER J1MIN,J2MAX      !calcula los limites reales del indice a medir
        INTEGER NCEFF
        INTEGER NCONTI,NABSOR
        REAL CA,CB,C3,C4
        REAL D1(NBDMAX),D2(NBDMAX),RL(NBDMAX),RG(NBDMAX),SUMRL
        REAL WLA,YDUM
        REAL TC,ETC                                              !flujo y error
        REAL FX(NBDMAX),EFX(NBDMAX)                !flujo y error en las bandas
        REAL SB,SR,ESB2,ESR2
        REAL SC(NCMAX),ESC2(NCMAX),MWB,MWR
        REAL S(NCMAX),ES(NCMAX),FSMEAN
        REAL FLUJOA,FLUJOB,XA,XB,XA_,XB_
        REAL COVAR_AA,COVAR_AB,COVAR_BB
        REAL SSIMUL(NCMAX)
        REAL WL(NCMAX),WL2(NCMAX)                          !pesos para el D4000
        REAL WVMIN,WVMAX           !limites reales en l.d.o. del indice a medir
        REAL F,FF
        REAL FRAC,PLAW(NCMAX)             !fraccion de la power law y power law
        REAL W0,WZ,MEAN1PLAW,MEAN2PLAW,KPLAW
        REAL SMAX,RINTERP
        REAL AMC,BMC       !parametros de recta en ajuste por minimos cuadrados
        REAL SUMX,SUMY,SUM0,SUMXY,SUMXX  !sumatorios para ajuste por min. cuad.
        REAL SIGMA2,DETER     !variable local para ajuste por minimos cuadrados
        DOUBLE PRECISION SMEAN
        LOGICAL IFCHAN(NCMAX)   !indica los canales usados para medir el indice
C
C PARAMETROS Y VARIABLES PARA LAS SIMULACIONES
        REAL PI2
        PARAMETER(PI2=6.283185307)
ccc     REAL MWC,CC,INTSC,INTB,INTC,INTR
        REAL SQ2ES(NCMAX),R1,R2,FINDEXSIMUL,SXX,SXXB,SXXR
        REAL WLA1,WLA2,COV
        REAL RANRED
C
        COMMON/BLKINDEX0/A,ERR
        COMMON/BLKINDEX1/FLUX
        COMMON/BLKINDEX2/WLMIN
        COMMON/BLKINDEX3/NSCAN,NCHAN,STWV,DISP
        COMMON/BLKINDEX4/NCREST
        COMMON/BLKINDEX8A/LPLOT,LPERRI           !dibujar o no, con/sin errores
        COMMON/BLKINDEX9/NSIMUL,NSEED
        COMMON/BLKINDEX10A/ALPHAPLAW,EBVPLAW
        COMMON/BLKINDEX10B/CBANDPLAW
        COMMON/BLKINDEX10C/RESCHAN
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C NOTA: el espectro a medir (scan numero I) se introduce en la matriz S(),
C       mientras que el error asociado (en caso necesario) se almacena en
C       la matriz ES(). Para todos los calculos se normaliza el valor de S()
C       (escalando igualmente ES()) para que los sumatorios introduzcan el
C       menor error posible.
C------------------------------------------------------------------------------
C Valores en caso de fallo
        FINDEX=0.
        EINDEX=0.
        EJJGG=0.
        SN=0.
        MIDEIND=1          !hasta que no se mida el indice, el retorno es fallo
C Evitamos WARNING de compilacion "may be used uninitialized in this function"
        ESB2=0.
        ESR2=0.
        SB=0.
        SR=0.
        MWB=0.
        MWR=0.
        DETER=0.
        SUM0=0.
        SUMX=0.
        SUMXX=0.
C------------------------------------------------------------------------------
C Protecciones
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: scan number out of range in subroutine '//
     +     'MIDEIND.'
          RETURN
        END IF
C
        IF((NCRES.LT.1).OR.(NCRES.GT.NCREST))THEN
          WRITE(*,101)'ERROR: invalid NCRES value in subroutine '//
     +     'MIDEIND.'
          RETURN
        END IF
C------------------------------------------------------------------------------
C Numero de bandas, limites y anchura teniendo en cuenta la velocidad radial
        NCONTI=0
        NABSOR=0
        IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
          NBAND=3             !indices atomicos y moleculares (lo mas probable)
        ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN
          NBAND=2                   !el D4000/B4000/color solo tiene dos bandas
        ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN !para los indices genericos:
          NCONTI=(ITI/100)                        !numero de bandas de continuo
          NABSOR=ITI-NCONTI*100                !numero de bandas de absorciones
          NBAND=NCONTI+NABSOR
        ELSEIF((ITI.LT.-1).AND.(ITI.GT.-100))THEN !indices pendiente
          NCONTI=-ITI
          NABSOR=0
          NBAND=-ITI
        ELSE
          WRITE(*,100) 'ITI: '
          WRITE(*,*) ITI
          WRITE(*,101) 'ERROR: invalid ITI value in subroutine MIDEIND.'
          RETURN
        END IF
C..............................................................................
        DO NB=1,NBAND                                          !para cada banda
          CA=WV(2*NB-1)*RCVEL1                           !redshifted wavelength
          CB=WV(2*NB)*RCVEL1                             !redshifted wavelength
          C3=(CA-WLMIN)/DISP+1.                           !band limit (channel)
          C4=(CB-WLMIN)/DISP                              !band limit (channel)
          IF((C3.LT.1.).OR.(C4.GT.REAL(NCHAN))) RETURN      !index out of range
          J1(NB)=INT(C3)                          !band limit: integer(channel)
          J2(NB)=INT(C4)                          !band limit: integer(channel)
          D1(NB)=C3-REAL(J1(NB))            !fraction (excess) of first channel
          D2(NB)=C4-REAL(J2(NB))             !fraction (defect) of last channel
          RL(NB)=CB-CA                       !redshifted band width (angstroms)
          RG(NB)=C4-C3                        !redshifted band width (channels)
        END DO
C..............................................................................
        J1MIN=J1(1)     !calculamos limites por si las bandas no estan en orden
        J2MAX=J2(1)
        WVMIN=WV(1)                                             !idem en l.d.o.
        WVMAX=WV(2)
        DO NB=2,NBAND
          IF(J1MIN.GT.J1(NB)) J1MIN=J1(NB)
          IF(J2MAX.LT.J2(NB)) J2MAX=J2(NB)
          IF(WVMIN.GT.WV(2*NB-1)) WVMIN=WV(2*NB-1)
          IF(WVMAX.LT.WV(2*NB)) WVMAX=WV(2*NB)
        END DO
C------------------------------------------------------------------------------
C Fijamos los canales a usar para medir el indice (utilizando la variable
C logica evitamos el problema de la posible superposicion de las bandas)
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.              !inicializamos: ningun canal utilizado
        END DO
C
        DO NB=1,NBAND                                    !recorremos las bandas
          DO J=J1(NB),J2(NB)+1
            IFCHAN(J)=.TRUE.
          END DO
        END DO
C
        NCEFF=0             !contamos el numero de canales effectivo a utilizar
        DO J=1,NCHAN
          IF(IFCHAN(J)) NCEFF=NCEFF+1
        END DO
C------------------------------------------------------------------------------
C Espectro (y errores) a medir (sumamos los scans que van desde NS1 a NS2)
        DO J=1,NCHAN
          S(J)=0.
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.
          END DO
        END IF
C
        DO I=NS1,NS2
          DO J=1,NCHAN
            S(J)=S(J)+A(J,I)
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO I=NS1,NS2
            DO J=1,NCHAN
              ES(J)=ES(J)+ERR(J,I)*ERR(J,I)
            END DO
          END DO
          DO J=1,NCHAN
            ES(J)=SQRT(ES(J))
          END DO
        END IF
C------------------------------------------------------------------------------
C Utilizacion de una curva respuesta no promedio: util para el calculo del
C error de calibracion en flujo (notar que la matriz A(J,I) se calibra en
C flujo al principio del programa y por eso hay que multiplicar por la curva
C respuesta inicial y dividir por la nueva curva respuesta).
        IF(NCRES.NE.1)THEN
          DO J=1,NCHAN
            S(J)=S(J)*FLUX(J,1)/FLUX(J,NCRES)
          END DO
C
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              ES(J)=ES(J)*FLUX(J,1)/FLUX(J,NCRES)
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C Si deseamos conocer el efecto de un continuo no termico de la forma:
C F(lambda) = k * (lambda)**(alpha-2.0), la variable FFPLAW sera un numero
C entero entre 0 y 10 (fraccion de luz entre 0.0 y 1.0). En caso contrario
C FFPLAW = -1 y no alteramos el espectro.
C Si introducimos una ley de potencias, el exponente es ALPHAPLAW. Hay que
C calcular la constante de la ley de potencias en funcion de la contribucion
C deseada. La fraccion de luz se refiere a la contribucion de la ley de
C potencias en los canales disponibles para medir la banda fotometrica elegida.
C Tenemos en cuenta el efecto del redshift en la ley de potencias y la
C extincion interna.
        IF((FFPLAW.GE.1).AND.(FFPLAW.LE.10))THEN
          DO J=1,NCHAN                          !calculamos la ley de potencias
            WZ=STWV+REAL(J-1)*DISP               !longitud de onda (redshift z)
            W0=WZ/RCVEL1                         !longitud de onda (redshift 0)
            PLAW(J)=
     +       (10.**(-.4*EBVPLAW*RINTERP(W0))*W0**(ALPHAPLAW-2.))/RCVEL1
          END DO
C..............................................................................
          IF(FFPLAW.EQ.10)THEN        !solo contribucion de la ley de potencias
            DO J=1,NCHAN
              S(J)=PLAW(J)
            END DO
            IF(CERR.EQ.'y')THEN
              DO J=1,NCHAN
                ES(J)=0.
              END DO
            END IF
C..............................................................................
C NOTA: en el caso siguiente el error no cambia (pero al aumentar S(J), 
C el error relativo disminuye, tendiendo a cero a medida que FRAC tiende a uno)
          ELSE                             !mezcla de objeto + ley de potencias
            SMAX=S(1)                                           !maximo de S(J)
            DO J=2,NCHAN
              IF(SMAX.LT.S(J)) SMAX=S(J)
            END DO
            MEAN1PLAW=0.     !numero de cuentas en la banda fotometrica elegida
            MEAN2PLAW=0.     !numero de cuentas en la banda fotometrica elegida
            DO J=1,NCHAN
              MEAN1PLAW=MEAN1PLAW+RESCHAN(J)*S(J)/SMAX
              MEAN2PLAW=MEAN2PLAW+RESCHAN(J)*PLAW(J)/SMAX
            END DO
            FRAC=REAL(FFPLAW)/10.              !fraccion de la ley de potencias
            KPLAW=(FRAC/(1.0-FRAC))*MEAN1PLAW/MEAN2PLAW  !cte. ley de potencias
            DO J=1,NCHAN
              S(J)=S(J)+KPLAW*PLAW(J)
            END DO
C..............................................................................
          END IF
        END IF
C------------------------------------------------------------------------------
C Normalizamos espectro (y errores) a medir, usando la senhal solo en la region
C del indice
        SMEAN=0.D0
        DO J=1,NCHAN
          IF(IFCHAN(J)) SMEAN=SMEAN+DBLE(S(J))      !solo canales en las bandas
        END DO
        SMEAN=SMEAN/DBLE(NCEFF)              !valor promedio (DOUBLE PRECISION)
        FSMEAN=REAL(SMEAN)                               !valor promedio (REAL)
C
        DO J=1,NCHAN
          S(J)=S(J)/FSMEAN                               !normalizamos espectro
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=ES(J)/FSMEAN                          !"normalizamos" errores
          END DO
        END IF
C------------------------------------------------------------------------------
C Senhal/ruido promedio en las bandas del indice
        IF((CERR.EQ.'y').AND.(FFPLAW.NE.10))THEN
          SN=0.                       !redundande, pero lo dejamos por claridad
          DO J=1,NCHAN
            IF(IFCHAN(J)) SN=SN+S(J)/ES(J)
          END DO
          SN=SN/REAL(NCEFF)
          SN=SN/SQRT(DISP)                  !damos la S/N promedio por Angstrom
        END IF
        IF(LONLY_SN) RETURN
C------------------------------------------------------------------------------
C Calculamos pseudo continuo en indices moleculares y atomicos (ITI=1,2)
C (formulas en Tesis de JJGG, pag. 35)
C
C NOTA: al transformar las integrales en sumatorios, habria que multiplicar
C cada valor de la funcion a integrar por el incremento en longitud de onda
C (que en el sumatorio coincide con DISP). Sin embargo, este valor es factor
C comun y puede salir fuera del sumatorio, por lo que solo hace falta
C introducirlo al final. Asimismo, como al calcular los limites de las bandas,
C J1() y J2(), hemos tenido presente la vel. radial, la anchura de las bandas,
C RL(), es mayor en un factor (1+z) que la anchura cuando el objeto no presenta
C velocidad radial. Es decir, el sumatorio se extiende sobre una region (en 
C longitud de onda) algo mayor. Sin embargo el efecto queda anulado al dividir
C por RL() que se encuentra ensanchado en el mismo factor.
C
        IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
C..............................................................................
          SB=0.                              !cuentas promedio en la banda azul
          DO J=J1(1),J2(1)+1                          !recorremos la banda azul
            IF(J.EQ.J1(1))THEN
              F=1.-D1(1)
            ELSEIF(J.EQ.J2(1)+1)THEN
              F=D2(1)
            ELSE
              F=1.
            END IF
            SB=SB+F*S(J)
          END DO
          SB=SB*DISP                                     !completamos sumatorio
          SB=SB/RL(1)                   !dividimos por anchura de la banda azul
C
          IF(CERR.EQ.'y')THEN
            ESB2=0.
            DO J=J1(1),J2(1)+1
              IF(J.EQ.J1(1))THEN
                F=1.-D1(1)
              ELSEIF(J.EQ.J2(1)+1)THEN
                F=D2(1)
              ELSE
                F=1.
              END IF
              ESB2=ESB2+F*F*ES(J)*ES(J)
            END DO
            ESB2=ESB2*DISP*DISP                          !completamos sumatorio
            ESB2=ESB2/(RL(1)*RL(1))     !dividimos por anchura de la banda azul
          END IF
C..............................................................................
          SR=0.                              !cuentas promedio en la banda roja
          DO J=J1(3),J2(3)+1                          !recorremos la banda roja
            IF(J.EQ.J1(3))THEN
              F=1.-D1(3)
            ELSEIF(J.EQ.J2(3)+1)THEN
              F=D2(3)
            ELSE
              F=1.
            END IF
            SR=SR+F*S(J)
          END DO
          SR=SR*DISP                                     !completamos sumatorio
          SR=SR/RL(3)                   !dividimos por anchura de la banda roja
C
          IF(CERR.EQ.'y')THEN
            ESR2=0.
            DO J=J1(3),J2(3)+1
              IF(J.EQ.J1(3))THEN
                F=1.-D1(3)
              ELSEIF(J.EQ.J2(3)+1)THEN
                F=D2(3)
              ELSE
                F=1.
              END IF
              ESR2=ESR2+F*F*ES(J)*ES(J)
            END DO
            ESR2=ESR2*DISP*DISP                          !completamos sumatorio
            ESR2=ESR2/(RL(3)*RL(3))     !dividimos por anchura de la banda roja
          END IF
C..............................................................................
C Trabajamos en la escala en l.d.o. sin corregir de Vrad (es decir, con las
C bandas de los indices desplazadas a Vrad correspondiente). Se obtiene lo
C mismo si mantenemos las bandas de los indices a Vrad=0 pero al calcular WLA
C para cada canal J dividimos por RCVEL1 para obtener la escala en l.d.o.
C corregida de Vrad.
          MWB=(WV(1)+WV(2))/2.             !mean wavelength blue band at Vrad=0
          MWB=MWB*RCVEL1                               !idem a Vrad considerado
          MWR=(WV(5)+WV(6))/2.              !mean wavelength red band at Vrad=0
          MWR=MWR*RCVEL1                               !idem a Vrad considerado
C Calculamos el valor del pseudo continuo desde el borde mas azul de todas las
C bandas al borde mas rojo de todas las bandas (no importa que las banda
C "azul" no sea la mas azul, etc.)
          DO J=J1MIN,J2MAX+1
            WLA=REAL(J-1)*DISP+STWV                         !l.d.o. del canal J
            SC(J)=(SB*(MWR-WLA)+SR*(WLA-MWB))/(MWR-MWB)
          END DO
C
          IF(CERR.EQ.'y')THEN
            DO J=J1MIN,J2MAX+1
              WLA=REAL(J-1)*DISP+STWV                       !l.d.o. del canal J
              ESC2(J)=ESB2*(MWR-WLA)*(MWR-WLA)+ESR2*(WLA-MWB)*(WLA-MWB)
              ESC2(J)=ESC2(J)/((MWR-MWB)*(MWR-MWB))
            END DO
          END IF
C..............................................................................
C En los dibujos escalamos todo a SMEAN
          IF(LPLOT)THEN                              !dibujamos pseudo continuo
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(7)
                CALL PGSLS(4)
              END DO
            END IF
            WLA=WVMIN*RCVEL1
            YDUM=SB*(MWR-WLA)/(MWR-MWB)+SR*(WLA-MWB)/(MWR-MWB)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(WLA,YDUM*FSMEAN)
            END DO
            WLA=WVMAX*RCVEL1
            YDUM=SB*(MWR-WLA)/(MWR-MWB)+SR*(WLA-MWB)/(MWR-MWB)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGDRAW(WLA,YDUM*FSMEAN)
            END DO
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
                CALL PGSLS(1)
              END DO
            END IF
          END IF
C..............................................................................
        END IF
C------------------------------------------------------------------------------
C Pesos para el D4000: debido a que la integral es  F(ldo)*d(nu), hay que
C multiplicar el flujo por el cuadrado de la longitud de onda, y de este
C modo la integral se transforma en F(ldo)*d(ldo).
C
C NOTA: estamos corrigiendo WLA de Vrad, aunque luego medimos el indice con
C los limites de las bandas multiplicados por (1+z). Esto no es importante
C porque solo hay un factor Cte=(1+z) en los pesos que no influye en el
C indice. Sin embargo, vamos a trabajar con los pesos asi porque de esta
C manera la normalizacion a 4000 siempre nos proporciona pesos cercanos a uno
C independientemente del valor de Vrad.
        IF(ITI.EQ.3)THEN
          DO J=1,NCHAN
            WLA=REAL(J-1)*DISP+STWV                         !l.d.o. del canal J 
            WLA=WLA/RCVEL1               !l.d.o. del canal J corrigiendo a Vrad
            WLA=WLA/4000.      !normalizamos a 4000 para tener todo proximo a 1
            WL(J)=WLA*WLA
          END DO
C
          IF(CERR.EQ.'y')THEN  !para los errores es mejor los pesos al cuadrado
            DO J=1,NCHAN
              WL2(J)=WL(J)*WL(J)
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C Para el B4000 usaremos codigo comun con el D4000, donde los pesos son iguales
C a uno.
        IF(ITI.EQ.4)THEN
          DO J=1,NCHAN
            WL(J)=1.0
          END DO
C
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              WL2(J)=1.0
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C Para el CO de Kleinmann and Hall (1986)
        IF(ITI.EQ.5)THEN
          DO J=1,NCHAN
            WL(J)=1.0
          END DO
C
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              WL2(J)=1.0
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C Para los indices genericos calculamos el pseudocontinuo ajustando por minimos
C cuadrados (pesando con errores si procede) a todos los pixels de las bandas 
C de continuo.
        IF((ITI.GE.101).AND.(ITI.LE.9999))THEN
C..............................................................................
C si no trabajamos con errores, hacemos todos estos iguales a uno para utilizar
C las mismas formulas
          IF(CERR.EQ.'n')THEN
            DO J=1,NCHAN
              ES(J)=1.0
            END DO
          END IF
C..............................................................................
C calculamos la recta del continuo mediante minimos cuadrados (la recta es
C de la forma y= amc * x + bmc
C (NOTA: para la variable "x" utilizamos el valor del numero de pixel en lugar
C de la longitud de onda porque, en principio, son numeros mas pequenhos)
          SUM0=0.
          SUMX=0.
          SUMY=0.
          SUMXY=0.
          SUMXX=0.
          DO NB=1,NCONTI               !recorremos todas las bandas de continuo
            DO J=J1(NB),J2(NB)+1  !recorremos todos los pixels de dichas bandas
              IF(J.EQ.J1(NB))THEN        !comprobamos efecto de borde izquierdo
                F=1.-D1(NB)
              ELSEIF(J.EQ.J2(NB)+1)THEN                !efecto de borde derecho
                F=D2(NB)
              ELSE                                  !pixels sin efecto de borde
                F=1.
              END IF
              SIGMA2=ES(J)*ES(J)
              SUM0=SUM0+F/SIGMA2
              SUMX=SUMX+F*REAL(J)/SIGMA2
              SUMY=SUMY+F*S(J)/SIGMA2
              SUMXY=SUMXY+F*REAL(J)*S(J)/SIGMA2
              SUMXX=SUMXX+F*REAL(J)*REAL(J)/SIGMA2
            END DO
          END DO
          DETER=SUM0*SUMXX-SUMX*SUMX
          AMC=(SUM0*SUMXY-SUMX*SUMY)/DETER
          BMC=(SUMXX*SUMY-SUMX*SUMXY)/DETER
C..............................................................................
C dibujamos el continuo calculado (escalando a SMEAN)
          IF(LPLOT)THEN
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(7)
                CALL PGSLS(4)
              END DO
            END IF
            YDUM=AMC*REAL(J1MIN)+BMC
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(STWV+REAL(J1MIN-1)*DISP,YDUM*FSMEAN)
            END DO
            YDUM=AMC*REAL(J2MAX)+BMC
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGDRAW(STWV+REAL(J2MAX-1)*DISP,YDUM*FSMEAN)
            END DO
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
                CALL PGSLS(1)
              END DO
            END IF
          END IF
C..............................................................................
C calculamos el pseudocontinuo desde el borde mas azul de todas las bandas
C hasta el borde mas rojo
          DO J=J1MIN,J2MAX+1
            SC(J)=AMC*REAL(J)+BMC
          END DO
C
          IF(CERR.EQ.'y')THEN
            DO J=J1MIN,J2MAX+1
              ESC2(J)=0.
              DO NB=1,NCONTI
                DO JJ=J1(NB),J2(NB)+1
                  SIGMA2=ES(JJ)*ES(JJ)
                  ESC2(J)=ESC2(J)+(
     +             (SUM0*REAL(JJ)/SIGMA2-SUMX/SIGMA2)*REAL(J)/DETER+
     +             (SUMXX/SIGMA2-SUMX*REAL(JJ)/SIGMA2)/DETER
     +             )**2 *SIGMA2
                END DO
              END DO
            END DO
          END IF
C..............................................................................
        END IF
C------------------------------------------------------------------------------
C Para los indices "pendiente" calculamos la recta ajustando por minimos
C cuadrados (pesando con errores si procede) a todos los pixels de las bandas 
C que definen el indice
        IF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
C..............................................................................
C si no trabajamos con errores, hacemos todos estos iguales a uno para utilizar
C las mismas formulas
          IF(CERR.EQ.'n')THEN
            DO J=1,NCHAN
              ES(J)=1.0
            END DO
          END IF
C..............................................................................
C calculamos la recta mediante minimos cuadrados (la recta es de la forma 
C y= amc * x + bmc )
C (NOTA: para la variable "x" utilizamos el valor del numero de pixel en lugar
C de la longitud de onda porque, en principio, son numeros mas pequenhos)
          SUM0=0.
          SUMX=0.
          SUMY=0.
          SUMXY=0.
          SUMXX=0.
          DO NB=1,NCONTI               !recorremos todas las bandas de continuo
            DO J=J1(NB),J2(NB)+1  !recorremos todos los pixels de dichas bandas
              IF(J.EQ.J1(NB))THEN        !comprobamos efecto de borde izquierdo
                F=1.-D1(NB)
              ELSEIF(J.EQ.J2(NB)+1)THEN                !efecto de borde derecho
                F=D2(NB)
              ELSE                                  !pixels sin efecto de borde
                F=1.
              END IF
              SIGMA2=ES(J)*ES(J)
              SUM0=SUM0+F/SIGMA2
              SUMX=SUMX+F*REAL(J)/SIGMA2
              SUMY=SUMY+F*S(J)/SIGMA2
              SUMXY=SUMXY+F*REAL(J)*S(J)/SIGMA2
              SUMXX=SUMXX+F*REAL(J)*REAL(J)/SIGMA2
            END DO
          END DO
          DETER=SUM0*SUMXX-SUMX*SUMX
          AMC=(SUM0*SUMXY-SUMX*SUMY)/DETER
          BMC=(SUMXX*SUMY-SUMX*SUMXY)/DETER
          IF(CERR.EQ.'y')THEN
            COVAR_AA=(SUM0*SUM0*SUMXX-SUM0*SUMX*SUMX)/DETER/DETER
            COVAR_AB=(SUMX*SUMX*SUMX-SUM0*SUMX*SUMXX)/DETER/DETER
            COVAR_BB=(SUMXX*SUMXX*SUM0-SUMXX*SUMX*SUMX)/DETER/DETER
          END IF
C..............................................................................
C dibujamos la recta calculada (escalando a SMEAN)
          IF(LPLOT)THEN
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(7)
                CALL PGSLS(4)
              END DO
            END IF
            YDUM=AMC*REAL(J1MIN)+BMC
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(STWV+REAL(J1MIN-1)*DISP,YDUM*FSMEAN)
            END DO
            YDUM=AMC*REAL(J2MAX)+BMC
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGDRAW(STWV+REAL(J2MAX-1)*DISP,YDUM*FSMEAN)
            END DO
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
                CALL PGSLS(1)
              END DO
            END IF
          END IF
C..............................................................................
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C medimos indices y errores
        IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN        !indices moleculares y atomicos
          TC=0.
          DO J=J1(2),J2(2)+1                       !recorremos la banda central
            IF(J.EQ.J1(2))THEN
              F=1.-D1(2)
            ELSEIF(J.EQ.J2(2)+1)THEN
              F=D2(2)
            ELSE
              F=1.
            END IF
            TC=TC+F*S(J)/SC(J)
          END DO
          TC=TC*DISP                                     !completamos sumatorio
          IF(ITI.EQ.1)THEN                                    !indice molecular
            FINDEX=-2.5*ALOG10(TC/RL(2))
          ELSEIF(ITI.EQ.2)THEN                                  !indice atomico
            FINDEX=RL(2)-TC
            FINDEX=FINDEX/RCVEL1                       !correccion por redshift
          END IF
C
          IF(CERR.EQ.'y')THEN             !los errores se suman cuadraticamente
            ETC=0.
            DO J=J1(2),J2(2)+1                     !recorremos la banda central
              IF(J.EQ.J1(2))THEN
                F=1.-D1(2)
              ELSEIF(J.EQ.J2(2)+1)THEN
                F=D2(2)
              ELSE
                F=1
              END IF
              ETC=ETC+F*F*(S(J)*S(J)*ESC2(J)+
     +         SC(J)*SC(J)*ES(J)*ES(J))/(SC(J)*SC(J)*SC(J)*SC(J))
C An~adimos terminos cruzados (covarianza)
              WLA1=REAL(J-1)*DISP+STWV
              DO JJ=J1(2),J2(2)+1
                IF(JJ.NE.J)THEN
                  IF(JJ.EQ.J1(2))THEN
                    FF=1.-D1(2)
                  ELSEIF(JJ.EQ.J2(2)+1)THEN
                    FF=D2(2)
                  ELSE
                    FF=1.
                  END IF
                  WLA2=REAL(JJ-1)*DISP+STWV
                  COV=(MWR-WLA1)*(MWR-WLA2)*ESB2+
     +             (WLA1-MWB)*(WLA2-MWB)*ESR2
                  COV=COV/((MWR-MWB)*(MWR-MWB))
ccc               type*,cov/sqrt(esc2(j)*esc2(jj))        !es muy proximo a uno
                  ETC=ETC+
     +             FF*F*S(J)*S(JJ)*COV/(SC(J)*SC(J)*SC(JJ)*SC(JJ))
                END IF
              END DO
            END DO
            ETC=SQRT(ETC)*DISP   !completamos sumatorio y hacemos raiz cuadrada
            IF(ITI.EQ.1)THEN                                  !indice molecular
              EINDEX=CTE/(10**(-.4*FINDEX))*ETC/RL(2)
            ELSEIF(ITI.EQ.2)THEN                                !indice atomico
              EINDEX=ETC
              EINDEX=EINDEX/RCVEL1                     !correccion por redshift
            END IF
          END IF
C------------------------------------------------------------------------------
        ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4).OR.(ITI.EQ.5))THEN   !D4000,B4000,color
C NOTA: los sumatorios TC y ETC no incluyen el factor DISP, que corresponderia
C con el incremento (diferencial en la integral), dado que al ser un factor
C constante tampoco altera el resultado a la hora de computar el D4000 (o el
C B4000, o el color).
          DO NB=1,NBAND                                        !bucle en bandas
            TC=0.
            DO J=J1(NB),J2(NB)+1                        !recorremos la banda NB
              IF(J.EQ.J1(NB))THEN
                F=1.-D1(NB)
              ELSEIF(J.EQ.J2(NB)+1)THEN
                F=D2(NB)
              ELSE
                F=1
              END IF
              TC=TC+F*S(J)*WL(J)
            END DO
            FX(NB)=TC
          END DO
          FINDEX=FX(2)/FX(1)
C
          IF(CERR.EQ.'y')THEN             !los errores se suman cuadraticamente
            DO NB=1,NBAND
              ETC=0.
              DO J=J1(NB),J2(NB)+1
                IF(J.EQ.J1(NB))THEN
                  F=1.-D1(NB)
                ELSEIF(J.EQ.J2(NB)+1)THEN
                  F=D2(NB)
                ELSE
                  F=1.
                END IF
                ETC=ETC+F*F*ES(J)*ES(J)*WL2(J)
              END DO
              EFX(NB)=ETC
            END DO
            EINDEX=FX(1)*FX(1)*EFX(2)+FX(2)*FX(2)*EFX(1)
            EINDEX=SQRT(EINDEX)/(FX(1)*FX(1))
          END IF
C
          IF(ITI.EQ.5)THEN
            FINDEX=FINDEX*RL(1)/RL(2)
            EINDEX=EINDEX*RL(1)/RL(2)
          END IF
C..............................................................................
C dibujamos en el D4000 (o el B4000) los valores medios en ambas bandas
          IF(LPLOT)THEN
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(7)
                CALL PGSLS(4)
              END DO
            END IF
C
            WLA=WV(1)*RCVEL1
            YDUM=0.
            DO J=J1(1),J2(1)+1
              YDUM=YDUM+S(J)
            END DO
            YDUM=YDUM/REAL(J2(1)-J1(1)+2)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(WLA,YDUM*FSMEAN)
            END DO
            WLA=WV(2)*RCVEL1
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGDRAW(WLA,YDUM*FSMEAN)
            END DO
C
            WLA=WV(3)*RCVEL1
            YDUM=0.
            DO J=J1(2),J2(2)+1
              YDUM=YDUM+S(J)
            END DO
            YDUM=YDUM/REAL(J2(2)-J1(2)+2)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMOVE(WLA,YDUM*FSMEAN)
            END DO
            WLA=WV(4)*RCVEL1
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGDRAW(WLA,YDUM*FSMEAN)
            END DO
C
            IF(NCRES.NE.1)THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(LCOLOR(ITERM)) CALL PGSCI(1)
                CALL PGSLS(1)
              END DO
            END IF
          END IF
C------------------------------------------------------------------------------
        ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
C recorremos las bandas con absorciones
          TC=0.
          SUMRL=0.
          DO NB=1,NABSOR
            DO J=J1(NCONTI+NB),J2(NCONTI+NB)+1     !recorremos todos los pixels
              IF(J.EQ.J1(NCONTI+NB))THEN !comprobamos efecto de borde izquierdo
                F=1.-D1(NCONTI+NB)
              ELSEIF(J.EQ.J2(NCONTI+NB)+1)THEN         !efecto de borde derecho
                F=D2(NCONTI+NB)
              ELSE                                  !pixels sin efecto de borde
                F=1.
              END IF
              TC=TC+F*FWV(NB)*S(J)/SC(J)    !multiplicamos absorcion por factor
            END DO
            SUMRL=SUMRL+FWV(NB)*RL(NCONTI+NB)
          END DO
          TC=TC*DISP                                  !completamos el sumatorio
          FINDEX=SUMRL-TC         !medimos el indice generico como los atomicos
          FINDEX=FINDEX/RCVEL1                         !correccion por redshift
C
          IF(CERR.EQ.'y')THEN             !los errores se suman cuadraticamente
            ETC=0.
            DO NB=1,NABSOR
              DO J=J1(NCONTI+NB),J2(NCONTI+NB)+1   !recorremos todos los pixels
                IF(J.EQ.J1(NCONTI+NB))THEN           !efecto de borde izquierdo
                  F=1.-D1(NCONTI+NB)
                ELSEIF(J.EQ.J2(NCONTI+NB)+1)THEN       !efecto de borde derecho
                  F=D2(NCONTI+NB)
                ELSE                                !pixels sin efecto de borde
                  F=1.
                END IF
                ETC=ETC+F*F*FWV(NB)*FWV(NB)*(S(J)*S(J)*ESC2(J)+
     +           SC(J)*SC(J)*ES(J)*ES(J))/(SC(J)*SC(J)*SC(J)*SC(J))
C An~adimos terminos cruzados (covarianza)
C Nota: los sumatorios ya fueron calculados mas arriba
                DO NBB=1,NABSOR
                  DO JJ=J1(NCONTI+NBB),J2(NCONTI+NBB)+1
                    IF((NBB.NE.NB).OR.(JJ.NE.J))THEN
                      IF(JJ.EQ.J1(NCONTI+NBB))THEN
                        FF=1.-D1(NCONTI+NBB)
                      ELSEIF(JJ.EQ.J2(NCONTI+NBB)+1)THEN
                        FF=D2(NCONTI+NBB)
                      ELSE
                        FF=1.
                      END IF
                      COV=REAL(J)*REAL(JJ)*
     +                 (SUM0*SUM0*SUMXX-SUM0*SUMX*SUMX)/(DETER*DETER)+
     +                 (REAL(J)+REAL(JJ))*
     +                 (SUMX*SUMX*SUMX-SUM0*SUMX*SUMXX)/(DETER*DETER)+
     +                 (SUMXX*SUMXX*SUM0-SUMXX*SUMX*SUMX)/(DETER*DETER)
                      ETC=ETC+FWV(NB)*FWV(NBB)*
     +                 FF*F*S(J)*S(JJ)*COV/(SC(J)*SC(J)*SC(JJ)*SC(JJ))
                    END IF
                  END DO
                END DO
              END DO
            END DO
            ETC=SQRT(ETC)*DISP   !completamos sumatorio y hacemos raiz cuadrada
            EINDEX=ETC/RCVEL1                           !correcion por redshift
          END IF
C------------------------------------------------------------------------------
        ELSEIF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
C calculamos el valor de la recta en las posicion central de la banda mas azul,
C XB, y de la banda mas roja, XA
          DO NB=1,NCONTI
            CA=WV(2*NB-1)*RCVEL1                         !redshifted wavelength
            CB=WV(2*NB)*RCVEL1                           !redshifted wavelength
            C3=(CA-WLMIN)/DISP+1.                         !band limit (channel)
            C4=(CB-WLMIN)/DISP                            !band limit (channel)
            IF(NB.EQ.1)THEN
              XA=(C3+C4)/2.0
              XB=(C3+C4)/2.0
            ELSE
              XA_=(C3+C4)/2.0
              XB_=(C3+C4)/2.0
              XA=AMAX1(XA,XA_)
              XB=AMIN1(XB,XB_)
            END IF
          END DO
          FLUJOA=AMC*XA+BMC
          FLUJOB=AMC*XB+BMC
          FINDEX=FLUJOA/FLUJOB
C
          IF(LPLOT)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(2)
              CALL PGPOINT(1,STWV+(XA-1.0)*DISP,FLUJOA*FSMEAN,17)
              IF(LCOLOR(ITERM)) CALL PGSCI(4)
              CALL PGPOINT(1,STWV+(XB-1.0)*DISP,FLUJOB*FSMEAN,17)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
C
          IF(CERR.EQ.'y')THEN
            EINDEX=SQRT(BMC*BMC*COVAR_AA+AMC*AMC*COVAR_BB-
     +       2.0*AMC*BMC*COVAR_AB)*(XA-XB)/((AMC*XB+BMC)*(AMC*XB+BMC))
          END IF
C
        END IF
C------------------------------------------------------------------------------
C Retorno con exito
        MIDEIND=0                       !el indice ha sido medido adecuadamente
ccc        WRITE(*,*)I,RCVEL1,NCRES,FINDEX,EINDEX
C------------------------------------------------------------------------------
C
C Con  la siguientes lineas comentadas saltamos el calculo de errores con 
C las formulas de JJGG's Thesis
C------------------------------------------------------------------------------
c
c calculo de errores usando las formulas de JJGG's Thesis
c        if((iti.eq.1).or.(iti.eq.2))then
c          if(cerr.eq.'y')then
c
c            mwc=(wv(3)+wv(4))/2.
c            mwc=mwc*rcvel1
c            cc=(sb*(mwr-mwc)+sr*(mwc-mwb))/(mwr-mwb)
c
c            intsc=0.
c            do j=j1(2),j2(2)
c              intsc=intsc+s(j)
c            end do
c            j=j1(2)
c            intsc=intsc-d1(2)*s(j)
c            j=j2(2)+1
c            intsc=intsc+d2(2)*s(j)
c            intsc=intsc*disp
c
c            intb=0.
c            do j=j1(1),j2(1)
c              intb=intb+s(j)*s(j)/(es(j)*es(j))
c            end do
c            j=j1(1)
c            intb=intb-d1(1)*s(j)*s(j)/(es(j)*es(j))
c            j=j2(1)+1
c            intb=intb+d2(1)*s(j)*s(j)/(es(j)*es(j))
ccc            intb=intb*disp !sobra
c            intb=sb*sb/intb
c
c            intc=0.
c            do j=j1(2),j2(2)
c              intc=intc+s(j)*s(j)/(es(j)*es(j))
c            end do
c            j=j1(2)
c            intc=intc-d1(2)*s(j)*s(j)/(es(j)*es(j))
c            j=j2(2)+1
c            intc=intc+d2(2)*s(j)*s(j)/(es(j)*es(j))
ccc            intc=intc*disp !sobra
c            intc=1./intc
c
c            intr=0.
c            do j=j1(3),j2(3)
c              intr=intr+s(j)*s(j)/(es(j)*es(j))
c            end do
c            j=j1(3)
c            intr=intr-d1(3)*s(j)*s(j)/(es(j)*es(j))
c            j=j2(3)+1
c            intr=intr+d2(3)*s(j)*s(j)/(es(j)*es(j))
ccc            intr=intr*disp !sobra
c            intr=sr*sr/intr
c
c            ejjgg=intc+intb/(cc*cc)*((mwr-mwc)/(mwr-mwb))**2+
c     +       intr/(cc*cc)*((mwb-mwc)/(mwr-mwb))**2
c            ejjgg=intsc*sqrt(ejjgg)/cc
c
c            if(iti.eq.1)then
c              ejjgg=cte/(10**(-.4*findex))*ejjgg/rl(2)
c            elseif(iti.eq.2)then
c              ejjgg=ejjgg/rcvel1                        !correccion de redshift
c            end if
c          end if
c        else        !JJGG no midio el D4000, i.e., no tiene formulas de errores
c          ejjgg=0.                             !redundante pero queda mas claro
c        end if
C
C------------------------------------------------------------------------------
C Calculo de errores usando simulaciones numericas
        IF(CERR.EQ.'y')THEN
C
          IF(NSIMUL.GT.1)THEN
            DO J=1,NCHAN
              SQ2ES(J)=1.414213562*ES(J)   !raiz de 2 por el error de la senhal
            END DO
c
c..............................................................................
C Para ITI=1 e ITI=2 podemos simular el error en el calculo del continuo de
C dos formas distintas. Podemos introducir ruido en cada pixel de las bandas
C laterales y calcular el valor promedio en la bandas azul y roja y, a partir de
C aqui calcular el continuo. Por otro lado, como conocemos los errores en los 
C valores promedio en las bandas azul y roja tambien podemos, directamente, 
C recalcular el continuo variando estos valores promedio, sin necesidad de 
C recalcularlos desde el principio. Ambos caminos producen el mismo resultado 
C (salvo errores estadisticos de las simulaciones).
            IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
              ESIMU=0.
              DO K=1,NSIMUL
ccc             SB=0.
ccc             DO J=J1(1),J2(1)+1
ccc               R1=RANRED(NSEED)
ccc               R2=RANRED(NSEED)
ccc               SXX=SQ2ES(J)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
ccc               IF(J.EQ.J1(1))THEN
ccc                 SB=SB+(1.-D1(1))*(S(J)+SXX)
ccc               ELSEIF(J.EQ.J2(1)+1)THEN
ccc                 SB=SB+D2(1)*(S(J)+SXX)
ccc               ELSE
ccc                 SB=SB+(S(J)+SXX)
ccc               END IF
ccc             END DO
ccc             SB=SB*DISP
ccc             SB=SB/RL(1)
c
ccc             SR=0.
ccc             DO J=J1(3),J2(3)+1
ccc               R1=RANRED(NSEED)
ccc               R2=RANRED(NSEED)
ccc               SXX=SQ2ES(J)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
ccc               IF(J.EQ.J1(3))THEN
ccc                 SR=SR+(1.-D1(3))*(S(J)+SXX)
ccc               ELSEIF(J.EQ.J2(3)+1)THEN
ccc                 SR=SR+D2(3)*(S(J)+SXX)
ccc               ELSE
ccc                 SR=SR+(S(J)+SXX)
ccc               END IF
ccc             END DO
ccc             SR=SR*DISP
ccc             SR=SR/RL(3)
c
ccc             DO J=J1(1),J2(3)+1
ccc               WLA=REAL(J-1)*DISP+STWV
ccc               SC(J)=(SB*(MWR-WLA)+SR*(WLA-MWB))/(MWR-MWB)
ccc             END DO
c
                R1=RANRED(NSEED)
                R2=RANRED(NSEED)
                SXXB=SQRT(2.*ESB2)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
                R1=RANRED(NSEED)
                R2=RANRED(NSEED)
                SXXR=SQRT(2.*ESR2)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
                DO J=J1(1),J2(3)+1
                  WLA=REAL(J-1)*DISP+STWV
                  SC(J)=((SB+SXXB)*(MWR-WLA)+(SR+SXXR)*(WLA-MWB))/
     +             (MWR-MWB)
                END DO
c
                TC=0.
                DO J=J1(2),J2(2)+1
                  R1=RANRED(NSEED)
                  R2=RANRED(NSEED)
                  SXX=SQ2ES(J)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
                  IF(J.EQ.J1(2))THEN
                    TC=TC+(1.-D1(2))*(S(J)+SXX)/SC(J)
                  ELSEIF(J.EQ.J2(2)+1)THEN
                    TC=TC+D2(2)*(S(J)+SXX)/SC(J)
                  ELSE
                    TC=TC+(S(J)+SXX)/SC(J)
                  END IF
                END DO
                TC=TC*DISP
                IF(ITI.EQ.1)THEN
                  FINDEXSIMUL=-2.5*ALOG10(TC/RL(2))
                ELSE!IF(ITI.EQ.2)THEN
                  FINDEXSIMUL=RL(2)-TC
                  FINDEXSIMUL=FINDEXSIMUL/RCVEL1       !correccion por redshift
                END IF
                ESIMU=ESIMU+(FINDEX-FINDEXSIMUL)*(FINDEX-FINDEXSIMUL)
              END DO
C
              IF(NSIMUL.GT.1)THEN
                ESIMU=SQRT(ESIMU/REAL(NSIMUL-1))
              ELSE
                ESIMU=0.
              END IF
C..............................................................................
            ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4))THEN
              ESIMU=0.
              DO K=1,NSIMUL
                DO NB=1,NBAND
                  TC=0.
                  DO J=J1(NB),J2(NB)+1
                    R1=RANRED(NSEED)
                    R2=RANRED(NSEED)
                    SXX=SQ2ES(J)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
                    IF(J.EQ.J1(NB))THEN
                      TC=TC+(1.-D1(NB))*(S(J)+SXX)*WL(J)
                    ELSEIF(J.EQ.J2(NB)+1)THEN
                      TC=TC+D2(NB)*(S(J)+SXX)*WL(J)
                    ELSE
                      TC=TC+(S(J)+SXX)*WL(J)
                    END IF
                  END DO
                  FX(NB)=TC
                END DO
                FINDEXSIMUL=FX(2)/FX(1)
                ESIMU=ESIMU+(FINDEX-FINDEXSIMUL)*(FINDEX-FINDEXSIMUL)
              END DO
              IF(NSIMUL.GT.1)THEN
                ESIMU=SQRT(ESIMU/REAL(NSIMUL-1))
              ELSE
                ESIMU=0.
              END IF
C..............................................................................
            ELSEIF(ITI.EQ.5)THEN
              ESIMU=0.0
              STOP 'Sorry: this option is not yet implemented'
C..............................................................................
            ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
              ESIMU=0.
C calculamos los sumatorios constantes
              SUM0=0.
              SUMX=0.
              SUMXX=0.
              DO NB=1,NCONTI
                DO J=J1(NB),J2(NB)+1
                  IF(J.EQ.J1(NB))THEN
                    F=1.-D1(NB)
                  ELSEIF(J.EQ.J2(NB)+1)THEN
                    F=D2(NB)
                  ELSE
                    F=1.
                  END IF
                  SIGMA2=ES(J)*ES(J)
                  SUM0=SUM0+F/SIGMA2
                  SUMX=SUMX+F*REAL(J)/SIGMA2
                  SUMXX=SUMXX+F*REAL(J)*REAL(J)/SIGMA2
                END DO
              END DO
              DETER=SUM0*SUMXX-SUMX*SUMX
c
              SUMRL=0.
              DO NB=1,NABSOR
                SUMRL=SUMRL+FWV(NB)*RL(NCONTI+NB)
              END DO
c
C iniciamos las diferentes simulaciones
              DO K=1,NSIMUL
                DO J=J1MIN,J2MAX+1
                  R1=RANRED(NSEED)
                  R2=RANRED(NSEED)
                  SXX=SQ2ES(J)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
                  SSIMUL(J)=S(J)+SXX
                END DO
                SUMY=0.
                SUMXY=0.
                DO NB=1,NCONTI
                  DO J=J1(NB),J2(NB)+1
                    IF(J.EQ.J1(NB))THEN
                      F=1.-D1(NB)
                    ELSEIF(J.EQ.J2(NB)+1)THEN
                      F=D2(NB)
                    ELSE
                      F=1.
                    END IF
                    SIGMA2=ES(J)*ES(J)
                    SUMY=SUMY+F*SSIMUL(J)/SIGMA2
                    SUMXY=SUMXY+F*REAL(J)*SSIMUL(J)/SIGMA2
                  END DO
                END DO
                AMC=(SUM0*SUMXY-SUMX*SUMY)/DETER
                BMC=(SUMXX*SUMY-SUMX*SUMXY)/DETER
                DO J=J1MIN,J2MAX+1
                  SC(J)=AMC*REAL(J)+BMC
                END DO
                TC=0.
                DO NB=1,NABSOR
                  DO J=J1(NCONTI+NB),J2(NCONTI+NB)+1
                    IF(J.EQ.J1(NCONTI+NB))THEN
                      F=1.-D1(NCONTI+NB)
                    ELSEIF(J.EQ.J2(NCONTI+NB)+1)THEN
                      F=D2(NCONTI+NB)
                    ELSE
                      F=1.
                    END IF
                    TC=TC+F*FWV(NB)*SSIMUL(J)/SC(J)
                  END DO
                END DO
                TC=TC*DISP
                FINDEXSIMUL=SUMRL-TC
                FINDEXSIMUL=FINDEXSIMUL/RCVEL1
                ESIMU=ESIMU+(FINDEX-FINDEXSIMUL)*(FINDEX-FINDEXSIMUL)
              END DO
C calculamos la desviacion tipica
              IF(NSIMUL.GT.1)THEN
                ESIMU=SQRT(ESIMU/REAL(NSIMUL-1))
              ELSE
                ESIMU=0.
              END IF
C..............................................................................
            ELSEIF((ITI.GT.-100).AND.(ITI.LT.-1))THEN
              ESIMU=0.
C calculamos los sumatorios constantes
              SUM0=0.
              SUMX=0.
              SUMXX=0.
              DO NB=1,NCONTI
                DO J=J1(NB),J2(NB)+1
                  IF(J.EQ.J1(NB))THEN
                    F=1.-D1(NB)
                  ELSEIF(J.EQ.J2(NB)+1)THEN
                    F=D2(NB)
                  ELSE
                    F=1.
                  END IF
                  SIGMA2=ES(J)*ES(J)
                  SUM0=SUM0+F/SIGMA2
                  SUMX=SUMX+F*REAL(J)/SIGMA2
                  SUMXX=SUMXX+F*REAL(J)*REAL(J)/SIGMA2
                END DO
              END DO
              DETER=SUM0*SUMXX-SUMX*SUMX
C iniciamos las diferentes simulaciones
              DO K=1,NSIMUL
                DO J=J1MIN,J2MAX+1
                  R1=RANRED(NSEED)
                  R2=RANRED(NSEED)
                  SXX=SQ2ES(J)*SQRT(-1.*ALOG(1.-R1))*COS(PI2*R2)
                  SSIMUL(J)=S(J)+SXX
                END DO
                SUMY=0.
                SUMXY=0.
                DO NB=1,NCONTI
                  DO J=J1(NB),J2(NB)+1
                    IF(J.EQ.J1(NB))THEN
                      F=1.-D1(NB)
                    ELSEIF(J.EQ.J2(NB)+1)THEN
                      F=D2(NB)
                    ELSE
                      F=1.
                    END IF
                    SIGMA2=ES(J)*ES(J)
                    SUMY=SUMY+F*SSIMUL(J)/SIGMA2
                    SUMXY=SUMXY+F*REAL(J)*SSIMUL(J)/SIGMA2
                  END DO
                END DO
                AMC=(SUM0*SUMXY-SUMX*SUMY)/DETER
                BMC=(SUMXX*SUMY-SUMX*SUMXY)/DETER
                FLUJOA=AMC*XA+BMC
                FLUJOB=AMC*XB+BMC
                FINDEXSIMUL=FLUJOA/FLUJOB
                ESIMU=ESIMU+(FINDEX-FINDEXSIMUL)*(FINDEX-FINDEXSIMUL)
              END DO
C calculamos la desviacion tipica
              IF(NSIMUL.GT.1)THEN
                ESIMU=SQRT(ESIMU/REAL(NSIMUL-1))
              ELSE
                ESIMU=0.
              END IF
            END IF
C..............................................................................
C numero de simulaciones insuficiente
          ELSE
            ESIMU=0.                           !redundante pero queda mas claro
          END IF
C
        END IF
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
C Fin de subrutina
        END
