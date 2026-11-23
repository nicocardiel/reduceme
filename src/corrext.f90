!------------------------------------------------------------------------------
! Version 24-November-2026                                    file: corrext.f90
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
! Program: corrext
! Classification: extinction correction
! Description: Corrects spectra from atmospheric and interstellar extinction.
!
!Comment
!
! Corrige de extincion atmosferica, conocida la masa de aire y el observatorio,
! y de extincion interestelar, conocido E(B-V). Realiza tambien el trabajo con
! los errores si se solicita.
! 
! Para compilar
! f77pgp corrext corrext.f Lred.a Lfutils.a
!
        PROGRAM CORREXT
        USE RED_M
        USE FSPEC_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER NPMAX
        PARAMETER(NPMAX=1000)
!
!
        INTEGER I,J,L,LRED
        INTEGER NPTOS
        REAL EXT(NCMAX)
        REAL, ALLOCATABLE :: A(:,:)
        REAL RLA,RINTERPLIN
        REAL F1,F2,F3,H
        REAL X(NPMAX),Y(NPMAX)
        REAL EBV
        REAL RV,K_LAMBDA
        CHARACTER*1 COK,CAT,CIN,CERR,CFUN,COPC
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,ERRFILE,OUTFILE
        CHARACTER*255 REDUCEMEDIR
        CHARACTER*255 EXTAFILE,EXTIFILE
        LOGICAL LOK
!
        COMMON/BLKINTER/X,Y,NPTOS
!------------------------------------------------------------------------------
        THISPROGRAM='corrext'
        CALL WELCOME('24-November-2026')
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        ! allocate with the actual image dimensions (static arrays
        ! were limited by NCMAX,NSMAX); initialised as before
        IF(ALLOCATED(A)) DEALLOCATE(A)
        ALLOCATE(A(NCHAN,NSCAN))
        A=0.0
        DO I=1,NSCAN
          READ(14) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(14)
        IF(CERR.EQ.'y')THEN    !abrimos fichero de errores pero no lo lemos aun
          WRITE(*,100,ADVANCE='NO')'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(13,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !...match
        END IF
!------------------------------------------------------------------------------
5       WRITE(CDUMMY,*)STWV
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'STWV = '//CDUMMY
        WRITE(CDUMMY,*)DISP
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'DISP = '//CDUMMY
        WRITE(*,100,ADVANCE='NO')'Are these values OK (y/n) '
        COK(1:1)=READC('y','yn')
        IF(COK.EQ.'n')THEN
          WRITE(*,100,ADVANCE='NO')'New STWV '
          WRITE(CDUMMY,*)STWV
          STWV=READF(CDUMMY)
          WRITE(*,100,ADVANCE='NO')'Nes DISP '
          WRITE(CDUMMY,*)DISP
          DISP=READF(CDUMMY)
          GOTO 5
        END IF
!
        DO I=1,NCHAN
          EXT(I)=1.0
        END DO
!
        CALL GET_ENVIRONMENT_VARIABLE('reduceme_dir',REDUCEMEDIR)
        LRED=TRUELEN(REDUCEMEDIR)
!
        WRITE(*,100,ADVANCE='NO')&
     &     'Correct from atmospheric extinction (y/n) '
        CAT(1:1)=READC('y','yn')
        IF(CAT.EQ.'y')THEN
          WRITE(*,101) '1 - Curve from file'
          WRITE(*,101) '2 - Semi-analytical function'
          WRITE(*,100,ADVANCE='NO') 'Option (1/2) '
          CFUN(1:1)=READC('2','12')
          IF(CFUN.EQ.'1')THEN
            WRITE(*,100,ADVANCE='NO')'Extinction curve file name '
            IF(LRED.GT.0)THEN
              EXTAFILE(1:LRED)=REDUCEMEDIR(1:LRED)
              EXTAFILE(LRED+1:)='/files/extlp.dat'
            ELSE
              EXTAFILE='extlp.dat'
            END IF
            EXTAFILE=INFILEX(15,EXTAFILE,0,0,.0,.0,3,.FALSE.)
            NPTOS=1
10          READ(15,*,END=12)X(NPTOS),Y(NPTOS)
            IF(NPTOS.GT.NPMAX) STOP 'FATAL ERROR: file too large.'
            NPTOS=NPTOS+1
            GOTO 10
12          NPTOS=NPTOS-1
            WRITE(*,110)'No. points read: ',NPTOS
            CLOSE(15)
            WRITE(*,100,ADVANCE='NO')'Mean air mass '
            WRITE(CDUMMY,*)AIRMASS
            AIRMASS=READF(CDUMMY)
            WRITE(*,100,ADVANCE='NO')'Interpolating function...'
            DO I=1,NCHAN
              RLA=STWV+REAL(I-1)*DISP
              EXT(I)=10.**(0.4*AIRMASS*RINTERPLIN(RLA))
            END DO
            WRITE(*,101)'   ...OK'
          ELSE
            WRITE(*,100,ADVANCE='NO') 'f1 '
            F1=READF('1.50')
            WRITE(*,100,ADVANCE='NO') 'f2 '
            F2=READF('1.25')
            WRITE(*,100,ADVANCE='NO') 'f3 '
            F3=READF('0.05')
            WRITE(*,100,ADVANCE='NO') 'H '
            H=READF('2.168')
            WRITE(*,100,ADVANCE='NO')'Mean air mass '
            WRITE(CDUMMY,*)AIRMASS
            AIRMASS=READF(CDUMMY)
            WRITE(*,100,ADVANCE='NO')'Interpolating function...'
            DO I=1,NCHAN
              RLA=(STWV+REAL(I-1)*DISP)/10000.0 !en micras
              EXT(I)=10.**(0.4*AIRMASS*FATMEXT(RLA,F1,F2,F3,H))
            END DO
            WRITE(*,101)'   ...OK'
          END IF
        END IF
!
        WRITE(*,100,ADVANCE='NO')&
     &     'Correct from interstellar extinction (y/n) '
        CIN(1:1)=READC('y','yn')
        IF(CIN.EQ.'y')THEN
          WRITE(*,101) '1 - Curve from file'
          WRITE(*,101) '2 - Pre-defined curve'
          WRITE(*,100,ADVANCE='NO') 'Option (1/2) '
          CFUN(1:1)=READC('2','12')
          IF(CFUN.EQ.'1')THEN
            WRITE(*,100,ADVANCE='NO')'Extinction curve file name '
            IF(LRED.GT.0)THEN
              EXTIFILE(1:LRED)=REDUCEMEDIR(1:LRED)
              EXTIFILE(LRED+1:)='/files/extint.dat'
            ELSE
              EXTIFILE='extint.dat'
            END IF
            EXTIFILE=INFILEX(16,EXTIFILE,0,0,.0,.0,3,.FALSE.)
            NPTOS=1
20          READ(16,*,END=22)X(NPTOS),Y(NPTOS)
            IF(NPTOS.GT.NPMAX) STOP 'FATAL ERROR: file too large.'
            NPTOS=NPTOS+1
            GOTO 20
22          NPTOS=NPTOS-1
            WRITE(*,110)'No. points read: ',NPTOS
            CLOSE(16)
            WRITE(*,100,ADVANCE='NO')'E(B-V)'
            EBV=READF('@')
            WRITE(*,100,ADVANCE='NO')'Interpolating function...'
            DO I=1,NCHAN
              RLA=STWV+REAL(I-1)*DISP
              EXT(I)=EXT(I)*(10.**(0.4*EBV*RINTERPLIN(RLA)))
            END DO
            WRITE(*,101)'   ...OK'
          ELSE
            WRITE(*,101) '(1) Galaxy: Savage & Mathis (1979)'
            WRITE(*,101) '(2) Galaxy: Seaton (1979)'
            WRITE(*,101) '(3) Galaxy: Cardelli, Clayton and Mathis'//&
     &       ' (1989) + O_Donnell (1994)'
            WRITE(*,101) '(4) Galaxy: Fitzpatrick (1999)'
            WRITE(*,101) '(a) LMC: Howarth (1983)'
            WRITE(*,101) '(b) LMC (30 Doradus): Fitzpatrick (1985)'
            WRITE(*,101) '(m) SMC: Prevot et al. (1984) + '//&
     &       'Bouchet et al. (1985)'
            WRITE(*,101) '(p) Starburst: Calzetti (1997)'
            WRITE(*,101) '(q) Starburst: Calzetti et al. (2000)'
            WRITE(*,100,ADVANCE='NO') 'Option (1/2/3/4/a/b/m/p/q) '
            COPC(1:1)=READC('4','1234abmpq')
            WRITE(*,100,ADVANCE='NO') 'R_V '
            RV=READF('3.1')
            WRITE(*,100,ADVANCE='NO')'E(B-V)'
            EBV=READF('@')
            WRITE(*,100,ADVANCE='NO')'Interpolating function...'
            DO I=1,NCHAN
              RLA=STWV+REAL(I-1)*DISP
              CALL SCOMPEXT(COPC,RLA,RV,K_LAMBDA,LOK)
              IF(.NOT.LOK)THEN
                WRITE(*,101) 'FATAL ERROR: LOK=.FALSE. in SCOMPEXT'
                STOP
              END IF
              EXT(I)=EXT(I)*(10.**(0.4*EBV*K_LAMBDA))
            END DO
            WRITE(*,101)'   ...OK'
          END IF
        END IF
!------------------------------------------------------------------------------
        DO I=1,NSCAN
          DO J=1,NCHAN
            A(J,I)=A(J,I)*EXT(J)
          END DO
        END DO
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(17,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(17) (A(J,I),J=1,NCHAN)
        END DO
        CALL CLOSEFILEX(17)
!------------------------------------------------------------------------------
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN                               !leemos fichero de errores
            READ(13) (A(J,I),J=1,NCHAN)
          END DO
          CLOSE(13)
          DO I=1,NSCAN
            DO J=1,NCHAN
              A(J,I)=A(J,I)*EXT(J)
            END DO
          END DO
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(18,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(18) (A(J,I),J=1,NCHAN)
          END DO
          CALL CLOSEFILEX(18)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
!
!******************************************************************************
!
        REAL FUNCTION RINTERPLIN(LAMBDA)
        IMPLICIT NONE
        INTEGER NPMAX
        PARAMETER(NPMAX=1000)
        REAL LAMBDA
        INTEGER I
        REAL X(NPMAX),Y(NPMAX)
        INTEGER NPTOS
        INTEGER N1,N2
        COMMON/BLKINTER/X,Y,NPTOS
!------------------------------------------------------------------------------
        N1=0
        N2=0
!
        DO I=1,NPTOS
          IF(X(I).LE.LAMBDA)N1=I
        END DO
        DO I=NPTOS,1,-1
          IF(X(I).GE.LAMBDA)N2=I
        END DO
        IF((N1.EQ.0).OR.(N2.EQ.0))THEN
          WRITE(*,100,ADVANCE='NO') 'LAMBDA='
          WRITE(*,*) LAMBDA
          STOP 'N1=0 or N2=0 in RINTERPLIN'
        ELSEIF(N1.EQ.N2)THEN
          RINTERPLIN=Y(N1)
        ELSE
          RINTERPLIN=Y(N1)+((LAMBDA-X(N1))/(X(N2)-X(N1)))*(Y(N2)-Y(N1))
        END IF
100     FORMAT(A)
        END
