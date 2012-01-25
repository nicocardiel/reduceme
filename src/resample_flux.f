C------------------------------------------------------------------------------
C Version 6-June-2003                                     file: resample_flux.f
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
C Program: resample_flux
C Classification: arithmetic & manipulations
C Description: Transforms an image with an initial wavelength calibration
C (linear or log-linear) into another image with a new and linear wavelength 
C calibration.
C
Comment
C
C Este programa transforma una imagen con una determinada STWV y DISP a otra
C imagen con una nueva STWV y DISP.
C
        PROGRAM RESAMPLE_FLUX
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J
        INTEGER NCHAN2
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL STWV2,DISP2
        REAL S(NCMAX),SS(NCMAX)
        REAL CRVAL,CDELT,CRPIX
        REAL WLMIN,WLMAX,WVMIN,WVMAX
        CHARACTER*1 CERR,COPC
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
C------------------------------------------------------------------------------
        THISPROGRAM='resample'
        CALL WELCOME('6-June-2003')
C imagen de entrada (y errores)
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(INFILE,ERRFILE)
          WRITE(*,100)'Error file name '
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C------------------------------------------------------------------------------
        WRITE(*,101) '1.- CTYPE=WAVE     (linear)'
        WRITE(*,101) '2.- CTYPE=WAVE-LOG (log10)'
        WRITE(*,101) '3.- CTYPE=WAVE-LOG (ln)'
        WRITE(*,101) '4.- CTYPE=wavenumber'
        WRITE(*,100) 'Option (1/2/3/4) '
        COPC(1:1)=READC('@','1234')
C------------------------------------------------------------------------------
C confirmamos la calibracion en l.d.o. del espectro de entrada
        IF((COPC.EQ.'1').OR.(COPC.EQ.'2').OR.(COPC.EQ.'3'))THEN
          WRITE(CDUMMY,*) STWV
          WRITE(*,100) 'STWV of input spectrum (Angstroms) '
          STWV=READF(CDUMMY)
          WRITE(CDUMMY,*) DISP
          WRITE(*,100) 'DISP of input spectrum (Angstroms) '
          DISP=READF(CDUMMY)
        ELSE
          WRITE(CDUMMY,*) STWV
          WRITE(*,100) 'STWV of input spectrum (wavenumber [cm^-1]) '
          STWV=READF(CDUMMY)
          STWV=STWV*1.E-8
          WRITE(CDUMMY,*) DISP
          WRITE(*,100) 'DISP of input spectrum (wavenumber [cm^-1]) '
          DISP=READF(CDUMMY)
          DISP=DISP*1.E-8
        END IF
C
        CRVAL=STWV
        CRPIX=1.0
        CDELT=DISP
C------------------------------------------------------------------------------
C nuevos parametros de la imagen de salida
        IF(COPC.EQ.'1')THEN
          STWV2=STWV
          DISP2=DISP
        ELSEIF(COPC.EQ.'2')THEN
          WLMIN=10.0**(CRVAL+CDELT*(0.5-CRPIX))
          WLMAX=10.0**(CRVAL+CDELT*(REAL(NCHAN)+0.5-CRPIX))
          DISP2=(WLMAX-WLMIN)/REAL(NCHAN)
          STWV2=WLMIN+0.5*DISP2
        ELSEIF(COPC.EQ.'3')THEN
          WLMIN=EXP(CRVAL+CDELT*(0.5-CRPIX))
          WLMAX=EXP(CRVAL+CDELT*(REAL(NCHAN)+0.5-CRPIX))
          DISP2=(WLMAX-WLMIN)/REAL(NCHAN)
          STWV2=WLMIN+0.5*DISP2
        ELSE
          WVMIN=(STWV-(0.5-CRPIX)*DISP)
          WVMAX=(STWV+REAL(NCHAN)*DISP+(0.5-CRPIX)*DISP)
          DISP2=((1./WVMIN-1./WVMAX))/REAL(NCHAN)
          STWV2=1.0/WVMAX+0.5*DISP2
        END IF
C
        WRITE(CDUMMY,*) STWV2
        WRITE(*,100)'New STWV in output spectrum (Angstroms) '
        STWV2=READF(CDUMMY)
C
        WRITE(CDUMMY,*) DISP2
        WRITE(*,100)'New DISP in output spectrum (Angstroms) '
        DISP2=READF(CDUMMY)
C
        WRITE(CDUMMY,*) NCHAN
        WRITE(*,100)'New NCHAN '
        NCHAN2=READILIM(CDUMMY,1,NCMAX)
C------------------------------------------------------------------------------
C creamos los nuevos espectros
        DO I=1,NSCAN
          DO J=1,NCHAN
            S(J)=A(J,I)
          END DO
          CALL ULOGREB(COPC,S,NCHAN,CRVAL,CRPIX,CDELT,
     +     SS,NCHAN2,STWV2,DISP2)
          DO J=1,NCHAN2
            A(J,I)=SS(J)
          END DO
        END DO
C
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              S(J)=ERR(J,I)
            END DO
            CALL ULOGREB(COPC,S,NCHAN,CRVAL,CRPIX,CDELT,
     +       SS,NCHAN2,STWV2,DISP2)
            DO J=1,NCHAN2
              ERR(J,I)=SS(J)
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
C salvamos la nueva imagen
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN2,STWV2,DISP2,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (A(J,I),J=1,NCHAN2)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          CALL GUESSEF(OUTFILE,ERRFILE)
          WRITE(*,100)'Output error file name '
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,NCHAN2,STWV2,DISP2,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (ERR(J,I),J=1,NCHAN2)
          END DO
          CLOSE(31)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
