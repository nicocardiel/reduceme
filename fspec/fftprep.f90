!------------------------------------------------------------------------------
! Version 20-November-2026                                    file: fttprep.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE FFTPREP(N,SP,NC1,NC2,LNORM,COSBELL,LFILT,KFILTER,LPLOT,CNAME)
!
! Input: N,SP,NC1,NC2,LNORM,COSBELL,LFILT,KFILTER,LPLOT,CNAME
! Output: SP
!
! This subroutine prepares spectrum SP for cross correlation. For this purpose
! the following steps are performed: 
! - normalization of the spectrum, using only the data in the range 
!   [SP(NC1),SP(NC2)]
! - extraction of the useful spectral region: [SP(NC1),SP(NC2)] is converted
!   into [SP(1),SP(NC2-NC1+1)]
! - multiplication of the spectrum by a cosine bell defined from 1 to NC2-NC1+1
! - zero padding from NC2-NC1+2 to N, where N=2**K, being K an integer
! - computation of the Fast Fourier Transform of SP(1),...,SP(N) and obtention
!   of the power spectrum
! - filtering of the power spectrum by applying the user defined filter KFILTER
! - computation of the inverse FFT of the filtered power spectrum and obtention
!   of the filtered spectrum SP(1),...,SP(N)
! - extraction of the initial spectral region: SP(1),...,SP(NC2-NC1+1)
! - multiplication by the cosine bell defined from 1 to NC2-NC1+1 (again!)
! - zero padding from NC2-NC1+2 to N (again!)
!
! INTEGER N -> number of points in output (must be a power of 2)
! REAL SP(N) -> spectrum to be prepared for cross correlation
! INTEGER NC1 -> index which indicates the first SP value to be employed
! INTEGER NC2 -> index which indicates the last SP value to be employed
! LOGICAL LNORM -> if .TRUE. SP is normalized prior any other manipulation
! REAL COSBELL(N) -> cosine bell (defined from point number 1 to NC2-NC1+1)
! LOGICAL LFILT -> if .TRUE. SP is filtered applying the filter KFILTER
! REAL KFILTER(N) -> filter (defined from point number 1 to N)
! LOGICAL LPLOT -> if .TRUE. plot intermediate plots
! CHARACTER*(*) CNAME -> string description for plots
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FFTPREP(N,SP,NC1,NC2,LNORM,COSBELL,LFILT,KFILTER,&
     &   LPLOT,CNAME)
        USE REDLIB_M
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
!
        INTEGER N
        REAL SP(N)
        INTEGER NC1,NC2
        LOGICAL LNORM
        REAL COSBELL(N)
        LOGICAL LFILT
        REAL KFILTER(N)
        LOGICAL LPLOT
        CHARACTER*(*) CNAME
!
        INTEGER NMAX
        PARAMETER (NMAX=8192)      !Si se cambia, cambiar tambien en subrutinas
!
        INTEGER J
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NWX,NWY
        REAL X(NMAX)
        REAL XR(NMAX),XI(NMAX)
        REAL FMODULO(NMAX),FASE(NMAX)
        REAL POWERLOG(NMAX)
        REAL FMEAN
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XV1,XV2,YV1,YV2,CH
        REAL XVT1,XVT2,YVT1,YVT2,DVX,DVY
        DOUBLE PRECISION DMEAN
        LOGICAL LCOLOR(MAX_ID_RED)
!
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
! una proteccion
        IF(N.GT.NMAX)THEN
          WRITE(*,100,ADVANCE='NO')'FATAL ERROR in subroutine SUBCORR: '
          WRITE(*,101)'N.GT.NMAX!'
          STOP
        END IF
!------------------------------------------------------------------------------
! si hay que dibujar, representamos los datos de entrada
        IF(LPLOT)THEN
          NWX=1
          NWY=4
          XVT1=0.
          XVT2=1.
          YVT1=0.05                   !dejamos un poco de hueco para PGIDEN_RED
          YVT2=1.00
          DVX=(XVT2-XVT1)/2.
          DVY=(YVT2-YVT1)/4.
          DO J=1,N
            X(J)=REAL(J)
          END DO
          XV1=XVT1+REAL(NWX-1)*DVX
          XV2=XV1+DVX
          YV1=YVT1+REAL(NWY-1)*DVY
          YV2=YV1+DVY
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGQCH(CH)
            CALL PGSCH(0.5)
            CALL PGPAGE
            CALL PGIDEN_RED
            CALL AUTOPLOT(N,X,SP,NC1,NC2,&
     &       'channel','no. counts',' ',&
     &       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &       0,.FALSE.,'BCNTS','BCNTS',&
     &       101,1,&
     &       XV1,XV2,YV1,YV2)
            CALL PGMTEXT('T',2.0,0.0,0.0,&
     &       'Original DATA in useful region')
            CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
          END DO
        END IF
!------------------------------------------------------------------------------
! si hay que normalizar el espectro, lo hacemos (calculamos el factor de
! normalizacion entre NC1 y NC2 exclusivamente)
        IF(LNORM)THEN
          DMEAN=0.D0
          DO J=NC1,NC2
            DMEAN=DMEAN+DBLE(SP(J))
          END DO
          DMEAN=DMEAN/DBLE(NC2-NC1+1)
          FMEAN=REAL(DMEAN)
          DO J=NC1,NC2
            SP(J)=SP(J)/FMEAN
          END DO
!
          IF(LPLOT)THEN
            NWX=NWX+1
            IF(NWX.EQ.3)THEN
              NWX=1
              NWY=NWY-1
            END IF
            XV1=XVT1+REAL(NWX-1)*DVX
            XV2=XV1+DVX
            YV1=YVT1+REAL(NWY-1)*DVY
            YV2=YV1+DVY
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL AUTOPLOT(N,X,SP,NC1,NC2,&
     &         'channel','no. counts',' ',&
     &         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &         0,.FALSE.,'BCNTS','BCNTS',&
     &         101,1,&
     &         XV1,XV2,YV1,YV2)
              CALL PGMTEXT('T',2.0,0.0,0.0,&
     &         'Normalized DATA in useful region')
              CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
            END DO
          END IF
        END IF
!------------------------------------------------------------------------------
! extraemos el espectro subconjunto a utilizar
        IF(NC1.NE.1)THEN
          DO J=1,NC2-NC1+1
            SP(J)=SP(J+NC1-1)
          END DO
        END IF
        IF(NC2-NC1+1.LT.N)THEN
          DO J=NC2-NC1+2,N
            SP(J)=0.
          END DO
        END IF
!
        IF(LPLOT)THEN
          NWX=NWX+1
          IF(NWX.EQ.3)THEN
            NWX=1
            NWY=NWY-1
          END IF
          XV1=XVT1+REAL(NWX-1)*DVX
          XV2=XV1+DVX
          YV1=YVT1+REAL(NWY-1)*DVY
          YV2=YV1+DVY
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL AUTOPLOT(N,X,SP,1,N,&
     &       'channel','no. counts',' ',&
     &       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &       0,.FALSE.,'BCNTS','BNTS',&
     &       101,1,&
     &       XV1,XV2,YV1,YV2)
            CALL PGSCI(5)
            CALL PGWINDOW(XMIN,XMAX,-0.05,1.05)
            CALL PGBOX(' ',0.0,0,'CMTS',0.0,0)
            CALL PGSLS(4)
            CALL PGBIN(NC2-NC1+1,X,COSBELL,.TRUE.)
            CALL PGSLS(1)
            CALL PGMTEXT('R',3.0,0.5,0.5,'cosine bell')
            CALL PGSCI(1)
            CALL PGMTEXT('T',2.0,0.0,0.0,&
     &       'Extending DATA with zero padding')
            CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
          END DO
        END IF
!------------------------------------------------------------------------------
! aplicamos la campana de coseno
        DO J=1,NC2-NC1+1
          SP(J)=SP(J)*COSBELL(J)
        END DO
!
        IF(LPLOT)THEN
          NWX=NWX+1
          IF(NWX.EQ.3)THEN
            NWX=1
            NWY=NWY-1
          END IF
          XV1=XVT1+REAL(NWX-1)*DVX
          XV2=XV1+DVX
          YV1=YVT1+REAL(NWY-1)*DVY
          YV2=YV1+DVY
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL AUTOPLOT(N,X,SP,1,N,&
     &       'channel','no. counts',' ',&
     &       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &       0,.FALSE.,'BCNTS','BCNTS',&
     &       101,1,&
     &       XV1,XV2,YV1,YV2)
            CALL PGMTEXT('T',2.0,0.0,0.0,&
     &       'Applying cosine bell in useful region')
            CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
          END DO
        END IF
!------------------------------------------------------------------------------
! si hay que filtrar, filtramos
        IF(LFILT)THEN
          DO J=1,NC2-NC1+1
            XR(J)=SP(J)                                    !parte real para FFT
            XI(J)=0.                                     !parte imaginaria nula
          END DO
          IF(N.GT.NC2-NC1+1)THEN                                  !zero padding
            DO J=NC2-NC1+2,N
              XR(J)=0.
              XI(J)=0.
            END DO
          END IF
          CALL CFFT(N,XR,XI,1) !............................................FFT
          DO J=1,N
            FMODULO(J)=SQRT(XR(J)*XR(J)+XI(J)*XI(J))
            FASE(J)=ATAN2(XI(J),XR(J))
          END DO
!
          IF(LPLOT)THEN
            DO J=1,N
              POWERLOG(J)=ALOG10(FMODULO(J))
            END DO
            NWX=NWX+1
            IF(NWX.EQ.3)THEN
              NWX=1
              NWY=NWY-1
            END IF
            XV1=XVT1+REAL(NWX-1)*DVX
            XV2=XV1+DVX
            YV1=YVT1+REAL(NWY-1)*DVY
            YV2=YV1+DVY
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL AUTOPLOT(N,X,POWERLOG,1,N,&
     &         'K axis: discrete transform domain','log\\d10\\u(P)',&
     &          'frequency (in units of the Nyquist frequency)',&
     &         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &         0,.FALSE.,'BNTS','BNTS',&
     &         101,1,&
     &         XV1,XV2,YV1,YV2)
              CALL PGBOX(' ',0.0,0,'C',0.0,0)
              CALL PGWINDOW(-0.10,2.10,YMIN,YMAX)
              CALL PGBOX('CMTS',0.0,0,' ',0.0,0)
              CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
              CALL PGMTEXT('T',2.0,0.0,0.0,'Power Spectrum')
              CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
              CALL PGSCI(5)
              CALL PGWINDOW(XMIN,XMAX,-0.05,1.05)
              CALL PGBOX(' ',0.0,0,'CMTS',0.0,0)
              CALL PGSLS(4)
              CALL PGBIN(N,X,KFILTER,.TRUE.)
              CALL PGSLS(1)
              CALL PGMTEXT('R',3.0,0.5,0.5,'filter')
              CALL PGSCI(1)
            END DO
          END IF
!
          DO J=1,N
            FMODULO(J)=FMODULO(J)*KFILTER(J)
          END DO
          DO J=1,N
            XR(J)=FMODULO(J)*COS(FASE(J))
            XI(J)=FMODULO(J)*SIN(FASE(J))
          END DO
          CALL CFFT(N,XR,XI,-1)
          DO J=1,NC2-NC1+1
            SP(J)=XR(J)
          END DO
!
          IF(LPLOT)THEN
            NWX=NWX+1
            IF(NWX.EQ.3)THEN
              NWX=1
              NWY=NWY-1
            END IF
            XV1=XVT1+REAL(NWX-1)*DVX
            XV2=XV1+DVX
            YV1=YVT1+REAL(NWY-1)*DVY
            YV2=YV1+DVY
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL AUTOPLOT(N,X,SP,1,NC2-NC1+1,&
     &         'channel','no. counts',' ',&
     &         .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &         0,.FALSE.,'BCNTS','BCNTS',&
     &         101,1,&
     &         XV1,XV2,YV1,YV2)
              CALL PGMTEXT('T',2.0,0.0,0.0,&
     &         'Filtered DATA in useful region')
              CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
            END DO
          END IF
        END IF
!------------------------------------------------------------------------------
! aplicamos otra vez la campana de coseno
        DO J=1,NC2-NC1+1
          SP(J)=SP(J)*COSBELL(J)
        END DO
!
        IF(LPLOT)THEN
          NWX=NWX+1
          IF(NWX.EQ.3)THEN
            NWX=1
            NWY=NWY-1
          END IF
          XV1=XVT1+REAL(NWX-1)*DVX
          XV2=XV1+DVX
          YV1=YVT1+REAL(NWY-1)*DVY
          YV2=YV1+DVY
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL AUTOPLOT(N,X,SP,1,NC2-NC1+1,&
     &       'channel','no. counts',' ',&
     &       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &       0,.FALSE.,'BCNTS','BCNTS',&
     &       101,1,&
     &       XV1,XV2,YV1,YV2)
            CALL PGMTEXT('T',2.0,0.0,0.0,&
     &       'Applying cosine bell in useful region')
            CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
          END DO
        END IF
!------------------------------------------------------------------------------
! zero padding otra vez
        IF(N.GT.NC2-NC1+1)THEN
          DO J=NC2-NC1+2,N
            SP(J)=0.
          END DO
        END IF
!
        IF(LPLOT)THEN
          NWX=NWX+1
          IF(NWX.EQ.3)THEN
            NWX=1
            NWY=NWY-1
          END IF
          XV1=XVT1+REAL(NWX-1)*DVX
          XV2=XV1+DVX
          YV1=YVT1+REAL(NWY-1)*DVY
          YV2=YV1+DVY
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL AUTOPLOT(N,X,SP,1,N,&
     &       'channel','no. counts',' ',&
     &       .TRUE.,XMIN,XMAX,YMIN,YMAX,0.05,&
     &       0,.FALSE.,'BCNTS','BCNTS',&
     &       101,1,&
     &       XV1,XV2,YV1,YV2)
            CALL PGMTEXT('T',2.0,0.0,0.0,'Applying zero padding')
            CALL PGMTEXT('T',2.0,1.0,1.0,CNAME)
            CALL PGSCH(CH)
          END DO
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
