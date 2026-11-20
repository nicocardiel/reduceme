!------------------------------------------------------------------------------
! Version 20-November-2026                                   file:adnscgrad.f90
! Adapted from adnsc by Javier Cenarro
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
! Program: adnscgrad
! Classification: arithmetic & manipulations
! Description: Add scans of an image, generating a new image with the same 
! NCHAN and variable NSCAN.
!
!Comment
!------------------------------------------------------------------------------
        PROGRAM ADNSCGRAD
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
!
        INTEGER I,J,L
        INTEGER N1,N2
        INTEGER NC1,NC2,NCEFF
        REAL A(NCMAX,NSMAX)
        REAL FSCALE,SCANEFF,SPROM(NSMAX),STOT,DISTARCSEC,CSCAN
        CHARACTER*1 CFHEAD2
        CHARACTER*15 CSCANBINNING
        CHARACTER*255 FILENAME,HEADFILE
        LOGICAL LOOP
        LOGICAL IFCHAN(NCMAX)
!------------------------------------------------------------------------------
        THISPROGRAM='adnscgrad'
        CALL WELCOME('20-November-2026')
        WRITE(*,*)
!
        WRITE(*,101) 'NOTE: This program simply computes effective '//&
     &   'radii for binned regions.'
        WRITE(*,101) '      You must use adnsc to generate the image'//&
     &   ' with the appropriate binning.'
        WRITE(*,*)
        WRITE(*,101) '      In addition, you must be careful when'//&
     &   ' computing binnings that include'
        WRITE(*,101) '      the central scan (this program assumes'//&
     &   ' r=0, instead of 1/4th of the'
        WRITE(*,101) '      projected pixel width).'
        WRITE(*,*)
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        FILENAME=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
! Pedimos datos relevantes
        WRITE(*,100,ADVANCE='NO') 'Galaxy center scan: '
        READ(*,*) CSCAN
        WRITE(*,100,ADVANCE='NO') 'Plate scale (arcsec/scan): '
        READ(*,*) FSCALE
! definimos los canales a utilizar para calcular el flujo promedio
5       WRITE(*,101) '* Define channel region to compute the '//&
     &   'average signal in each spectrum:'
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE.
        END DO
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(*,100,ADVANCE='NO') 'Channel region NC1,NC2 (0,0=EXIT) '
          CALL READ2I('0,0',NC1,NC2)
          IF((NC1.EQ.0).AND.(NC2.EQ.0))THEN
            LOOP=.FALSE.
          ELSE
            IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
              WRITE(*,101) 'ERROR: invalid entry. Try again.'
            ELSE
              DO J=NC1,NC2
                IFCHAN(J)=.TRUE.
              END DO
            END IF
          END IF
        END DO
        !contamos los canales fuera del bucle anterior por si acaso el
        !usuario duplica o solapa regiones por equivocacion
        NCEFF=0
        DO J=1,NCHAN
          IF(IFCHAN(J))THEN
            NCEFF=NCEFF+1
          END IF
        END DO
        IF(NCEFF.EQ.0)THEN
          WRITE(*,101) 'ERROR: you must define a channel region!'
          WRITE(*,100,ADVANCE='NO') 'Press <CR> to continue...'
          READ(*,*)
          GOTO 5
        END IF
! Decidimos si salvamos los resultados en un fichero externos
        WRITE(*,100,ADVANCE='NO')&
     &     'Create log file (with binning regions) '//&
     &   'for index.f (y/n) '
        CFHEAD2(1:1)=READC('n','yn')
        IF(CFHEAD2.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Log file name'
          HEADFILE=OUTFILEX(27,'@',0,0,0.,0.,3,.FALSE.)
        END IF
!------------------------------------------------------------------------------
10      WRITE(*,100,ADVANCE='NO')'Scan region (0,0=EXIT) '
        CALL READ2I('0,0',N1,N2)
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 20
        IF((N1.LT.0).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
          WRITE(*,101)'ERROR: invalid entry. Try again.'
          GOTO 10
        END IF
!
        SCANEFF=0.
        STOT=0.
        DO I=N1,N2
           SPROM(I)=0.
           DO J=1,NCHAN
             IF(IFCHAN(J))THEN
               SPROM(I)=SPROM(I)+A(J,I) !flujo promedio en scan i
             END IF
           END DO
           SPROM(I)=SPROM(I)/REAL(NCEFF) !flujo promedio en regiones definidas
           SCANEFF=SCANEFF+SPROM(I)*REAL(I) !pesamos scan I con peso SPROM
           STOT=STOT+SPROM(I) !flujo total en el bin de scans (suma de pesos)
        END DO
        SCANEFF=SCANEFF/STOT
        DISTARCSEC=(SCANEFF-CSCAN)*FSCALE
!
        WRITE(*,100,ADVANCE='NO') '>>>'
        WRITE(*,100,ADVANCE='NO') '>>> Binning from scan '
        write(*,'(I5,A9,I5)') N1,' to scan ',N2
        WRITE(*,100,ADVANCE='NO')& 
     &   '>>> Effective (luminosity-weighted) scan..........: '
        WRITE(*,'(F9.2)') SCANEFF
        WRITE(*,100,ADVANCE='NO')& 
     &   '>>> Effective distance from galaxy center (arcsec): '
        WRITE(*,'(F9.2)') DISTARCSEC
!
        IF(CFHEAD2.EQ.'y')THEN
          WRITE(CSCANBINNING,'(I7.0,A1,I7.0)')N1,',',N2
          CALL RMBLANK(CSCANBINNING,CSCANBINNING,L)
          WRITE(27,'(A15,2(1X,F9.2))') CSCANBINNING(1:L),&
     &     SCANEFF,DISTARCSEC
        END IF
        GOTO 10
!------------------------------------------------------------------------------
20      IF(CFHEAD2.EQ.'y') CLOSE(27)
!
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
