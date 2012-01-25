C------------------------------------------------------------------------------
C Version 06-April-2001                                         file: leeposs.f
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
C Program: leeposs
C Classification: miscellany
C Description: Reads FITS images created with getimage from the Digital Sky
C Survey, allowing the determination of accurate coordinates and plotting of
C marks.
C
Comment
C
C Permite leer un fichero en formato FITS generado por POSS, dibuja la
C imagen y calcula la posicion (A.R. y DEC.) de los pixels marcados con
C la ayuda de un raton. Permite hacer un binning en la imagen.
C
        PROGRAM LEEPOSS
C
        IMPLICIT NONE
C
        INTEGER NSMAX,NCMAX
C ATENCION! Si se cambian los parametros hay que hacerlo tambien en las
C subrutinas.
        PARAMETER (NSMAX=2400,NCMAX=2400)
C
        INTEGER NSCAN,NCHAN
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
C
        INTEGER IHEADER
        INTEGER I,J,L
        INTEGER I1
        INTEGER J1,II,JJ
        INTEGER NPIXELS,NROWS
        INTEGER BITPIX,NAXIS1,NAXIS2
        INTEGER CNPIX1,CNPIX2
        INTEGER IEXTRAER
        INTEGER PLTRAH,PLTRAM
        INTEGER PLTDECD,PLTDECM
        INTEGER NBIN
        INTEGER NTERM,IDN(8)
        INTEGER IUNIT,BLOCKSIZE,ISTATUS,NKEYS,NSPACE,NAXIS(0:2)
        INTEGER FIRSTPIX,JROW(NCMAX),NFOUND,IREADWRITE
        REAL A(NCMAX,NSMAX)
        REAL SUM
        REAL XV1,XV2,YV1,YV2
        DOUBLE PRECISION PLTRAS,PLTDECS
        DOUBLE PRECISION PPO3,PPO6
        DOUBLE PRECISION XPIXELSZ,YPIXELSZ
        DOUBLE PRECISION AMDX1,AMDX2,AMDX3,AMDX4,AMDX5,AMDX6,AMDX7
        DOUBLE PRECISION AMDX8,AMDX9,AMDX10,AMDX11,AMDX12,AMDX13
        DOUBLE PRECISION AMDY1,AMDY2,AMDY3,AMDY4,AMDY5,AMDY6,AMDY7
        DOUBLE PRECISION AMDY8,AMDY9,AMDY10,AMDY11,AMDY12,AMDY13
        DOUBLE PRECISION DEXTRAER
        CHARACTER*1 CEXTRAER,CBIN
        CHARACTER*1 PLTDECSN
        CHARACTER*50 COMMENT
        CHARACTER*80 FILENAME,CLINEA
        CHARACTER*2880 BLOQUE
        LOGICAL LOGFILE
        LOGICAL LBITPIX,LNAXIS1,LNAXIS2
        LOGICAL LPLTRAH,LPLTRAM,LPLTRAS
        LOGICAL LPLTDECSN,LPLTDECD,LPLTDECM,LPLTDECS
        LOGICAL LPPO3,LPPO6
        LOGICAL LCOLOR(8)
        LOGICAL LXPIXELSZ,LYPIXELSZ
        LOGICAL LAMDX1,LAMDX2,LAMDX3,LAMDX4,LAMDX5,LAMDX6,LAMDX7
        LOGICAL LAMDX8,LAMDX9,LAMDX10,LAMDX11,LAMDX12,LAMDX13
        LOGICAL LAMDY1,LAMDY2,LAMDY3,LAMDY4,LAMDY5,LAMDY6,LAMDY7
        LOGICAL LAMDY8,LAMDY9,LAMDY10,LAMDY11,LAMDY12,LAMDY13
        LOGICAL LCNPIX1,LCNPIX2
        LOGICAL LROW(NCMAX),ANYNULL
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLKNBIN/NBIN
        COMMON/BLKPOSS1/PPO3,PPO6,XPIXELSZ,YPIXELSZ
        COMMON/BLKPOSS2/AMDX1,AMDX2,AMDX3,AMDX4,AMDX5,AMDX6,AMDX7,
     +   AMDX8,AMDX9,AMDX10,AMDX11,AMDX12,AMDX13
        COMMON/BLKPOSS3/AMDY1,AMDY2,AMDY3,AMDY4,AMDY5,AMDY6,AMDY7,
     +   AMDY8,AMDY9,AMDY10,AMDY11,AMDY12,AMDY13
        COMMON/BLKPOSS4/CNPIX1,CNPIX2
        COMMON/BLKPOSS5/PLTRAH,PLTRAM,PLTDECD,PLTDECM
        COMMON/BLKPOSS6/PLTRAS,PLTDECS
        COMMON/BLKPOSS7/PLTDECSN
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
10      WRITE(*,100)'Input file name'
        FILENAME(1:80)=READC('@','@')
        INQUIRE(FILE=FILENAME,EXIST=LOGFILE)
        IF(LOGFILE)THEN
          OPEN(14,FILE=FILENAME,STATUS='OLD',FORM='FORMATTED',
ccc     +     RECORDTYPE='FIXED',RECL=2880,ACCESS='DIRECT',ERR=999)
     +     RECL=2880,ACCESS='DIRECT',ERR=999)
        ELSE
          WRITE(*,101)'ERROR: This file does not exist. Try again.'
          GOTO 10
        END IF
C looking for header
        IHEADER=1
        LBITPIX=.FALSE.
        LNAXIS1=.FALSE.
        LNAXIS2=.FALSE.
        LPLTRAH=.FALSE.
        LPLTRAM=.FALSE.
        LPLTRAS=.FALSE.
        LPLTDECSN=.FALSE.
        LPLTDECD=.FALSE.
        LPLTDECM=.FALSE.
        LPLTDECS=.FALSE.
        LPPO3=.FALSE.
        LPPO6=.FALSE.
        LXPIXELSZ=.FALSE.
        LYPIXELSZ=.FALSE.
        LAMDX1=.FALSE.
        LAMDX2=.FALSE.
        LAMDX3=.FALSE.
        LAMDX4=.FALSE.
        LAMDX5=.FALSE.
        LAMDX6=.FALSE.
        LAMDX7=.FALSE.
        LAMDX8=.FALSE.
        LAMDX9=.FALSE.
        LAMDX10=.FALSE.
        LAMDX11=.FALSE.
        LAMDX12=.FALSE.
        LAMDX13=.FALSE.
        LAMDY1=.FALSE.
        LAMDY2=.FALSE.
        LAMDY3=.FALSE.
        LAMDY4=.FALSE.
        LAMDY5=.FALSE.
        LAMDY6=.FALSE.
        LAMDY7=.FALSE.
        LAMDY8=.FALSE.
        LAMDY9=.FALSE.
        LAMDY10=.FALSE.
        LAMDY11=.FALSE.
        LAMDY12=.FALSE.
        LAMDY13=.FALSE.
        LCNPIX1=.FALSE.
        LCNPIX2=.FALSE.
20      READ(14,'(A2880)',REC=IHEADER)BLOQUE
C looking for bytes/pixel and image dimensions
        IF(INDEX(BLOQUE,'BITPIX').NE.0)THEN
          LBITPIX=.TRUE.
          BITPIX=IEXTRAER(BLOQUE,'BITPIX')
          WRITE(*,100)'BITPIX = '
          WRITE(*,*)BITPIX
        END IF
        IF(INDEX(BLOQUE,'NAXIS1 ').NE.0)THEN
          LNAXIS1=.TRUE.
          NAXIS1=IEXTRAER(BLOQUE,'NAXIS1 ')
          WRITE(*,100)'NAXIS1 = '
          WRITE(*,*)NAXIS1
        END IF
        IF(INDEX(BLOQUE,'NAXIS2 ').NE.0)THEN
          LNAXIS2=.TRUE.
          NAXIS2=IEXTRAER(BLOQUE,'NAXIS2 ')
          WRITE(*,100)'NAXIS2 = '
          WRITE(*,*)NAXIS2
        END IF
C lookin for POSS parameters
        IF(INDEX(BLOQUE,'PLTRAH').NE.0)THEN
          LPLTRAH=.TRUE.
          PLTRAH=IEXTRAER(BLOQUE,'PLTRAH')
        END IF
        IF(INDEX(BLOQUE,'PLTRAM').NE.0)THEN
          LPLTRAM=.TRUE.
          PLTRAM=IEXTRAER(BLOQUE,'PLTRAM')
        END IF
        IF(INDEX(BLOQUE,'PLTRAS').NE.0)THEN
          LPLTRAS=.TRUE.
          PLTRAS=DEXTRAER(BLOQUE,'PLTRAS')
        END IF
        IF(INDEX(BLOQUE,'PLTDECSN').NE.0)THEN
          LPLTDECSN=.TRUE.
          PLTDECSN=CEXTRAER(BLOQUE,'PLTDECSN')
        END IF
        IF(INDEX(BLOQUE,'PLTDECD').NE.0)THEN
          LPLTDECD=.TRUE.
          PLTDECD=IEXTRAER(BLOQUE,'PLTDECD')
        END IF
        IF(INDEX(BLOQUE,'PLTDECM').NE.0)THEN
          LPLTDECM=.TRUE.
          PLTDECM=IEXTRAER(BLOQUE,'PLTDECM')
        END IF
        IF(INDEX(BLOQUE,'PLTDECS ').NE.0)THEN
          LPLTDECS=.TRUE.
          PLTDECS=DEXTRAER(BLOQUE,'PLTDECS ')
        END IF
        IF(INDEX(BLOQUE,'PPO3').NE.0)THEN
          LPPO3=.TRUE.
          PPO3=DEXTRAER(BLOQUE,'PPO3')
        END IF
        IF(INDEX(BLOQUE,'PPO6').NE.0)THEN
          LPPO6=.TRUE.
          PPO6=DEXTRAER(BLOQUE,'PPO6')
        END IF
        IF(INDEX(BLOQUE,'XPIXELSZ').NE.0)THEN
          LXPIXELSZ=.TRUE.
          XPIXELSZ=DEXTRAER(BLOQUE,'XPIXELSZ')
        END IF
        IF(INDEX(BLOQUE,'YPIXELSZ').NE.0)THEN
          LYPIXELSZ=.TRUE.
          YPIXELSZ=DEXTRAER(BLOQUE,'YPIXELSZ')
        END IF
        IF(INDEX(BLOQUE,'AMDX1 ').NE.0)THEN
          LAMDX1=.TRUE.
          AMDX1=DEXTRAER(BLOQUE,'AMDX1 ')
        END IF
        IF(INDEX(BLOQUE,'AMDX2 ').NE.0)THEN
          LAMDX2=.TRUE.
          AMDX2=DEXTRAER(BLOQUE,'AMDX2 ')
        END IF
        IF(INDEX(BLOQUE,'AMDX3 ').NE.0)THEN
          LAMDX3=.TRUE.
          AMDX3=DEXTRAER(BLOQUE,'AMDX3 ')
        END IF
        IF(INDEX(BLOQUE,'AMDX4').NE.0)THEN
          LAMDX4=.TRUE.
          AMDX4=DEXTRAER(BLOQUE,'AMDX4')
        END IF
        IF(INDEX(BLOQUE,'AMDX5').NE.0)THEN
          LAMDX5=.TRUE.
          AMDX5=DEXTRAER(BLOQUE,'AMDX5')
        END IF
        IF(INDEX(BLOQUE,'AMDX6').NE.0)THEN
          LAMDX6=.TRUE.
          AMDX6=DEXTRAER(BLOQUE,'AMDX6')
        END IF
        IF(INDEX(BLOQUE,'AMDX7').NE.0)THEN
          LAMDX7=.TRUE.
          AMDX7=DEXTRAER(BLOQUE,'AMDX7')
        END IF
        IF(INDEX(BLOQUE,'AMDX8').NE.0)THEN
          LAMDX8=.TRUE.
          AMDX8=DEXTRAER(BLOQUE,'AMDX8')
        END IF
        IF(INDEX(BLOQUE,'AMDX9').NE.0)THEN
          LAMDX9=.TRUE.
          AMDX9=DEXTRAER(BLOQUE,'AMDX9')
        END IF
        IF(INDEX(BLOQUE,'AMDX10').NE.0)THEN
          LAMDX10=.TRUE.
          AMDX10=DEXTRAER(BLOQUE,'AMDX10')
        END IF
        IF(INDEX(BLOQUE,'AMDX11').NE.0)THEN
          LAMDX11=.TRUE.
          AMDX11=DEXTRAER(BLOQUE,'AMDX11')
        END IF
        IF(INDEX(BLOQUE,'AMDX12').NE.0)THEN
          LAMDX12=.TRUE.
          AMDX12=DEXTRAER(BLOQUE,'AMDX12')
        END IF
        IF(INDEX(BLOQUE,'AMDX13').NE.0)THEN
          LAMDX13=.TRUE.
          AMDX13=DEXTRAER(BLOQUE,'AMDX13')
        END IF
        IF(INDEX(BLOQUE,'AMDY1 ').NE.0)THEN
          LAMDY1=.TRUE.
          AMDY1=DEXTRAER(BLOQUE,'AMDY1 ')
        END IF
        IF(INDEX(BLOQUE,'AMDY2 ').NE.0)THEN
          LAMDY2=.TRUE.
          AMDY2=DEXTRAER(BLOQUE,'AMDY2 ')
        END IF
        IF(INDEX(BLOQUE,'AMDY3 ').NE.0)THEN
          LAMDY3=.TRUE.
          AMDY3=DEXTRAER(BLOQUE,'AMDY3 ')
        END IF
        IF(INDEX(BLOQUE,'AMDY4').NE.0)THEN
          LAMDY4=.TRUE.
          AMDY4=DEXTRAER(BLOQUE,'AMDY4')
        END IF
        IF(INDEX(BLOQUE,'AMDY5').NE.0)THEN
          LAMDY5=.TRUE.
          AMDY5=DEXTRAER(BLOQUE,'AMDY5')
        END IF
        IF(INDEX(BLOQUE,'AMDY6').NE.0)THEN
          LAMDY6=.TRUE.
          AMDY6=DEXTRAER(BLOQUE,'AMDY6')
        END IF
        IF(INDEX(BLOQUE,'AMDY7').NE.0)THEN
          LAMDY7=.TRUE.
          AMDY7=DEXTRAER(BLOQUE,'AMDY7')
        END IF
        IF(INDEX(BLOQUE,'AMDY8').NE.0)THEN
          LAMDY8=.TRUE.
          AMDY8=DEXTRAER(BLOQUE,'AMDY8')
        END IF
        IF(INDEX(BLOQUE,'AMDY9').NE.0)THEN
          LAMDY9=.TRUE.
          AMDY9=DEXTRAER(BLOQUE,'AMDY9')
        END IF
        IF(INDEX(BLOQUE,'AMDY10').NE.0)THEN
          LAMDY10=.TRUE.
          AMDY10=DEXTRAER(BLOQUE,'AMDY10')
        END IF
        IF(INDEX(BLOQUE,'AMDY11').NE.0)THEN
          LAMDY11=.TRUE.
          AMDY11=DEXTRAER(BLOQUE,'AMDY11')
        END IF
        IF(INDEX(BLOQUE,'AMDY12').NE.0)THEN
          LAMDY12=.TRUE.
          AMDY12=DEXTRAER(BLOQUE,'AMDY12')
        END IF
        IF(INDEX(BLOQUE,'AMDY13').NE.0)THEN
          LAMDY13=.TRUE.
          AMDY13=DEXTRAER(BLOQUE,'AMDY13')
        END IF
        IF(INDEX(BLOQUE,'CNPIX1').NE.0)THEN
          LCNPIX1=.TRUE.
          CNPIX1=IEXTRAER(BLOQUE,'CNPIX1')
        END IF
        IF(INDEX(BLOQUE,'CNPIX2').NE.0)THEN
          LCNPIX2=.TRUE.
          CNPIX2=IEXTRAER(BLOQUE,'CNPIX2')
        END IF
C
C looking for end of header
        IF(INDEX(BLOQUE,' END    ').NE.0) GOTO 30
        IHEADER=IHEADER+1
        GOTO 20
C
30      WRITE(*,100)'Number of header blocks: '
        WRITE(*,*)IHEADER
C
        IF(LBITPIX)THEN
          IF(BITPIX.NE.16)THEN
            WRITE(*,101)'FATAL ERROR #2: BITPIX.NE.16'
            STOP
          END IF
        ELSE
          WRITE(*,101)'BITPIX not found'
        END IF
C
        IF((LNAXIS1).AND.(LNAXIS2))THEN
          NPIXELS=NAXIS1
          NROWS=NAXIS2
        ELSE
          WRITE(*,101)'WARNING: NAXIS1 and NAXIS2 have not been found'
          WRITE(*,*)
          WRITE(*,100)'Number of channels/frame'
          NPIXELS=READI('@')
          WRITE(*,100)'Number of scans/frame   '
          NROWS=READI('@')
        END IF
        IF(NPIXELS.GT.NCMAX)THEN
          WRITE(*,101)'FATAL ERROR #3: number of channels too large.'
          STOP
        END IF
        IF(NROWS.GT.NSMAX)THEN
          WRITE(*,101)'FATAL ERROR #4: number of scans too large.'
          STOP
        END IF
        NCHAN=NPIXELS
        NSCAN=NROWS
C
        IF(LCNPIX1.AND.LCNPIX2)THEN
C          WRITE(*,101)'CNPIX1 and CNPIX2 found'
        ELSE
          WRITE(*,101)'FATAL ERROR: POSS parameters not found.'
          STOP
        END IF
C
        IF(LPPO3.AND.LPPO6)THEN
C          WRITE(*,101)'PPO3 and PPO4 found'
        ELSE
          WRITE(*,101)'FATAL ERROR: POSS parameters not found.'
          STOP
        END IF
C
        IF(LXPIXELSZ.AND.LYPIXELSZ)THEN
C          WRITE(*,101)'XPIXELSZ and YPIXELSZ found'
        ELSE
          WRITE(*,101)'FATAL ERROR: POSS parameters not found.'
          STOP
        END IF
C
        IF(LAMDX1.AND.LAMDX2.AND.LAMDX3.AND.LAMDX4.AND.LAMDX5.AND.
     +   LAMDX6.AND.LAMDX7.AND.LAMDX8.AND.LAMDX9.AND.LAMDX10.AND.
     +   LAMDX11.AND.LAMDX12.AND.LAMDX13)THEN
C          WRITE(*,101)'AMDX1...AMDX13 found'
        ELSE
          WRITE(*,101)'FATAL ERROR: POSS parameters not found.'
          STOP
        END IF
C
        IF(LAMDY1.AND.LAMDY2.AND.LAMDY3.AND.LAMDY4.AND.LAMDY5.AND.
     +   LAMDY6.AND.LAMDY7.AND.LAMDY8.AND.LAMDY9.AND.LAMDY10.AND.
     +   LAMDY11.AND.LAMDY12.AND.LAMDY13)THEN
C          WRITE(*,101)'AMDY1...AMDY13 found'
        ELSE
          WRITE(*,101)'FATAL ERROR: POSS parameters not found.'
          STOP
        END IF
C
        IF(LPLTRAH.AND.LPLTRAM.AND.LPLTRAS.AND.LPLTDECSN.AND.
     +   LPLTDECD.AND.LPLTDECM.AND.LPLTDECS)THEN
C          WRITE(*,101)'R.A. and DEC. (plate center) found'
        ELSE
          WRITE(*,101)'FATAL ERROR: POSS parameters not found.'
          STOP
        END IF
C
        WRITE(*,101)'All the POSS parameters have been found.'
C
        CLOSE(14) !cerramos fichero para ahora volver a abrirlo y leer imagen
C------------------------------------------------------------------------------
        IUNIT=14
        IREADWRITE=0
        ISTATUS=0
C abrimos el fichero
        CALL FTOPEN(IUNIT,FILENAME,IREADWRITE,BLOCKSIZE,ISTATUS)
        CALL PRINTERROR(ISTATUS)
C determinamos el numero de keywords en la cabecera y las mostramos
        CALL FTGHSP(IUNIT,NKEYS,NSPACE,ISTATUS)
        print*,'ftghsp: ',iunit,nkeys,nspace,istatus
        CALL PRINTERROR(ISTATUS)
        DO I=1,NKEYS
          CALL FTGREC(IUNIT,I,CLINEA,ISTATUS)
          CALL PRINTERROR(ISTATUS)
          L=TRUELEN(CLINEA)
          WRITE(*,101)CLINEA(1:L)
        END DO
        IF(ISTATUS.EQ.0)THEN                                  !todo ha ido bien
          WRITE(*,101)'END'
          WRITE(*,*)
        END IF
C leemos BITPIX
        CALL FTGKYJ(IUNIT,'BITPIX',BITPIX,COMMENT,ISTATUS)
        CALL PRINTERROR(ISTATUS)
C comprobamos que NAXIS=2
        CALL FTGKYJ(IUNIT,'NAXIS',NAXIS(0),COMMENT,ISTATUS)
        CALL PRINTERROR(ISTATUS)
        IF(NAXIS(0).NE.2)THEN
          WRITE(*,101)'FATAL ERROR: NAXIS is not equal to 2.'
          WRITE(*,*) NAXIS(0)
          CALL FTCLOS(IUNIT,ISTATUS)
          CALL PRINTERROR(ISTATUS)
          STOP
        END IF
C leemos NAXIS1 y NAXIS2 [notar que el quinto parametro es NAXIS(1) en lugar
C de NAXIS para asi recuperar NAXIS(1) y NAXIS(2)]
        CALL FTGKNJ(IUNIT,'NAXIS',1,2,NAXIS(1),NFOUND,ISTATUS)
        CALL PRINTERROR(ISTATUS)
        IF(NAXIS(1).GT.NCMAX)THEN
          WRITE(*,101)'* FATAL ERROR in subroutine LEEFITS:'
          WRITE(*,101)'NAXIS(1) > NCMAX'
          STOP
        END IF
        IF(NAXIS(2).GT.NSMAX)THEN
          WRITE(*,101)'* FATAL ERROR in subroutine LEEFITS:'
          WRITE(*,101)'NAXIS(2) > NSMAX'
          STOP
        END IF
C leemos la imagen
        IF(BITPIX.EQ.16)THEN
          DO I=1,NAXIS(2)
            FIRSTPIX=(I-1)*NAXIS(1)+1
            CALL FTGPFJ(IUNIT,1,FIRSTPIX,NAXIS(1),JROW(1),LROW(1),
     +       ANYNULL,ISTATUS)
            CALL PRINTERROR(ISTATUS)
            DO J=1,NAXIS(1)
              A(J,I)=REAL(JROW(J))
            END DO
          END DO
        ELSE
          CLOSE(IUNIT)
          STOP
        END IF
C------------------------------------------------------------------------------
ccc40      WRITE(*,101)'OK! File read.'
        WRITE(*,101)'OK! File read.'
        CLOSE(IUNIT)
C
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL BUTTSYB(3)
        CALL BUTTQPR(XV1,XV2,YV1,YV2)
        CALL BUTTSPR(XV1,XV2,YV1,0.75)
        CALL BUTTQBR(XV1,XV2,YV1,YV2)
        CALL BUTTSBR(XV1,XV2,0.84,YV2)
        CALL BUTTSIT(.TRUE.)
        CALL BUTTSCH(0.8)
        WRITE(*,100)'Binning (y/n) '
        CBIN(1:1)=READC('n','yn')
        IF(CBIN.EQ.'y')THEN
          WRITE(*,100)'Bin size (1=no binning) '
          NBIN=READILIM('1',1,9999)
          NCHAN=NCHAN/NBIN
          NSCAN=NSCAN/NBIN
          DO I=1,NSCAN
            DO J=1,NCHAN
              I1=(I-1)*NBIN+1
              J1=(J-1)*NBIN+1
              SUM=0.
              DO II=I1,I1+NBIN-1
                DO JJ=J1,J1+NBIN-1
                  SUM=SUM+A(JJ,II)
                END DO
              END DO
              A(J,I)=SUM/REAL(NBIN*NBIN)
            END DO
          END DO
        ELSE
          NBIN=1
        END IF
        CALL SUBLOOK(FILENAME)
C
999     WRITE(*,101)'FATAL ERROR#0: I/O ERROR OPENING THE FILE.'
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C****************************************************************************
C
        INTEGER FUNCTION IEXTRAER(CADENA,SUBCADENA)
C Esta funcion retorna el numero entero que existe en la cadena CADENA
C despues de encontrar SUBCADENA y un simbolo igual (=)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
        CHARACTER*(*) SUBCADENA
C
        INTEGER POS0,POS1
        INTEGER L,I
C
        L=LEN(CADENA)
        POS0=INDEX(CADENA,SUBCADENA)
        IF((POS0.EQ.0).OR.(POS0.EQ.L))THEN
          IEXTRAER=0
          WRITE(*,101)'IEXTRAER WARNING #1'
          RETURN
        END IF
C looking for '='
        DO I=POS0+1,L
          IF(CADENA(I:I).EQ.'=') GOTO 10
        END DO
        IEXTRAER=0
        WRITE(*,101)'IEXTRAER WARNING #2'
        RETURN
C looking for a non-blank
10      POS0=I
        DO I=POS0+1,LEN(CADENA)
          IF(CADENA(I:I).NE.' ') GOTO 20
        END DO
        IEXTRAER=0
        WRITE(*,101)'IEXTRAER WARNING #3'
        RETURN
C looking for the end of the number
20      POS0=I
        DO I=POS0+1,LEN(CADENA)
          IF(CADENA(I:I).EQ.' ') GOTO 30
        END DO
        POS1=LEN(CADENA)
        READ(CADENA(POS0:POS1),*)IEXTRAER
        WRITE(*,101)'IEXTRAER WARNING #4'
        RETURN
30      POS1=I-1
        READ(CADENA(POS0:POS1),*)IEXTRAER
        RETURN
101     FORMAT(A)
        END
C
C****************************************************************************
C
        CHARACTER*1 FUNCTION CEXTRAER(CADENA,SUBCADENA)
C Esta funcion retorna el caracter que existe en la cadena CADENA
C despues de encontrar SUBCADENA y un simbolo igual (=)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
        CHARACTER*(*) SUBCADENA
C
        INTEGER POS0,POS1
        INTEGER L,I
C
        L=LEN(CADENA)
        POS0=INDEX(CADENA,SUBCADENA)
        IF((POS0.EQ.0).OR.(POS0.EQ.L))THEN
          CEXTRAER=' '
          WRITE(*,101)'CEXTRAER WARNING #1'
          RETURN
        END IF
C looking for '='
        DO I=POS0+1,L
          IF(CADENA(I:I).EQ.'=') GOTO 10
        END DO
        CEXTRAER=' '
        WRITE(*,101)'CEXTRAER WARNING #2'
        RETURN
C looking for a non-blank
10      POS0=I
        DO I=POS0+1,LEN(CADENA)
          IF(CADENA(I:I).NE.' ') GOTO 20
        END DO
        CEXTRAER=' '
        WRITE(*,101)'CEXTRAER WARNING #3'
        RETURN
C looking for the end of the string
20      POS0=I
        DO I=POS0+1,LEN(CADENA)
          IF(CADENA(I:I).EQ.' ') GOTO 30
        END DO
        POS1=LEN(CADENA)
        CEXTRAER=CADENA(POS0+1:POS0+1)
        WRITE(*,101)'CEXTRAER WARNING #4'
        RETURN
30      POS1=I-1
        CEXTRAER=CADENA(POS0+1:POS0+1)
        RETURN
101     FORMAT(A)
        END
C
C****************************************************************************
C
        DOUBLE PRECISION FUNCTION DEXTRAER(CADENA,SUBCADENA)
C Esta funcion retorna el numero real que existe en la cadena CADENA
C despues de encontrar SUBCADENA y un simbolo igual (=)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
        CHARACTER*(*) SUBCADENA
C
        INTEGER POS0,POS1
        INTEGER L,I
C
        L=LEN(CADENA)
        POS0=INDEX(CADENA,SUBCADENA)
        IF((POS0.EQ.0).OR.(POS0.EQ.L))THEN
          DEXTRAER=0.D0
          WRITE(*,101)'DEXTRAER WARNING #1'
          RETURN
        END IF
C looking for '='
        DO I=POS0+1,L
          IF(CADENA(I:I).EQ.'=') GOTO 10
        END DO
        DEXTRAER=0.D0
        WRITE(*,101)'DEXTRAER WARNING #2'
        RETURN
C looking for a non-blank
10      POS0=I
        DO I=POS0+1,LEN(CADENA)
          IF(CADENA(I:I).NE.' ') GOTO 20
        END DO
        DEXTRAER=0.D0
        WRITE(*,101)'DEXTRAER WARNING #3'
        RETURN
C looking for the end of the number
20      POS0=I
        DO I=POS0+1,LEN(CADENA)
          IF(CADENA(I:I).EQ.' ') GOTO 30
        END DO
        POS1=LEN(CADENA)
        READ(CADENA(POS0:POS1),*)DEXTRAER
        WRITE(*,101)'DEXTRAER WARNING #4'
        RETURN
30      POS1=I-1
        READ(CADENA(POS0:POS1),*)DEXTRAER
        RETURN
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C
C Permite dibujar la imagen completa o hacer zoom sobre ella
C Los parametros que se pasan a traves del COMMON son
C A(j,i) - la matriz imagen con i:scans, j:canales
C NSCAN: el numero de scans de la imagen
C NCHAN: el numero de canales de la imagen
        SUBROUTINE SUBLOOK(FILENAME)
        IMPLICIT NONE
        CHARACTER*(*) FILENAME
C
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INTEGER NSMAX,NCMAX
        PARAMETER (NSMAX=2400,NCMAX=2400)
        INTEGER NSCAN,NCHAN
ccc        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        REAL READF
C
        INTEGER NMAXOBJ                    !cambiar tambien en otras subrutinas
        PARAMETER (NMAXOBJ=8000)
C
        INTEGER NOBJ,N0
        INTEGER NC1,NC2,NS1,NS2
        INTEGER II,JJ
        INTEGER IIMIN,IIMAX,JJMIN,JJMAX
        INTEGER IXC1,IXC2,IYC1,IYC2
        INTEGER NTPT
        INTEGER NB,NBLOCAL
        INTEGER NBIN
        INTEGER NTERM,IDN(8),ITERM
        INTEGER IDOLD,IDNEW,PGOPEN
        REAL A(NCMAX,NSMAX)
        REAL BG,FG,FDUM
        REAL TR(6),XC,YC
        REAL MAXVAL,MINVAL
        INTEGER ARH,ARM
        REAL ARS
        INTEGER DECD,DECM
        REAL DECS
        REAL EQUINOXFILE
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XMINA,XMAXA,YMINA,YMAXA
        DOUBLE PRECISION MEANVAL,SIGMA
        CHARACTER*1 DECSIG,CQUIT,CIOMARKS,CMERGE,CJ2000
        CHARACTER*1 CH,CNEXT
        CHARACTER*50 CDUMMY
        CHARACTER*75 MARKFILE,OUTFILE
        LOGICAL INSIDE1
        LOGICAL LCOLOR(8),LOGFILE,LBEXIST
C Variables para poner las marcas
        INTEGER NTYPE(NMAXOBJ)
        REAL XOBJ(NMAXOBJ),YOBJ(NMAXOBJ),ROBJ(NMAXOBJ)               !en pixels
        REAL XOBJ_1,XOBJ_2,YOBJ_1,YOBJ_2
        REAL PAOBJ(NMAXOBJ),PA_12                                !PA in degrees
        REAL SWARCOBJ(NMAXOBJ),SWOBJ(NMAXOBJ)      !slit width (arcsec, pixels)
        INTEGER ARHOBJ(NMAXOBJ),ARHOBJ_1,ARHOBJ_2
        INTEGER ARMOBJ(NMAXOBJ),ARMOBJ_1,ARMOBJ_2
        REAL ARSOBJ(NMAXOBJ),ARSOBJ_1,ARSOBJ_2
        CHARACTER*1 DECSIGOBJ(NMAXOBJ),DECSIGOBJ_1,DECSIGOBJ_2
        INTEGER DECDOBJ(NMAXOBJ),DECDOBJ_1,DECDOBJ_2
        INTEGER DECMOBJ(NMAXOBJ),DECMOBJ_1,DECMOBJ_2
        REAL DECSOBJ(NMAXOBJ),DECSOBJ_1,DECSOBJ_2
        REAL RARCOBJ(NMAXOBJ)          !tamanho de la marca en segundos de arco
        REAL LASTRADIUS,LASTWIDTH
        CHARACTER*20 TEXTOBJ(NMAXOBJ)
        REAL WINSIZE,FNC1,FNC2,FNS1,FNS2
        CHARACTER*1 CWSIZE,CSELMARK,CRADIUS
        INTEGER NCOLTYPE(10)                    !colores asignados a las marcas
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLKNBIN/NBIN
        COMMON/BLKSYMB1/NOBJ
        COMMON/BLKSYMB2/XOBJ,YOBJ,ROBJ,PAOBJ,SWOBJ
        COMMON/BLKSYMB3/NTYPE
        COMMON/BLKSYMB4/ARHOBJ,ARMOBJ,ARSOBJ
        COMMON/BLKSYMB5/DECSIGOBJ,DECDOBJ,DECMOBJ,DECSOBJ
        COMMON/BLKSYMB6/RARCOBJ,SWARCOBJ
        COMMON/BLKSYMB7/TEXTOBJ
        COMMON/BLKSYMB8/NC1,NC2,NS1,NS2
        COMMON/BLKSYMB9/NCOLTYPE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        NOBJ=0                      !numero de objetos marcados sobre la imagen
        CJ2000='y'
        NCOLTYPE(1)=2                           !colores asignados a las marcas
        NCOLTYPE(2)=4
        NCOLTYPE(3)=5
        NCOLTYPE(4)=7
        NCOLTYPE(5)=7
        NCOLTYPE(6)=6
        NCOLTYPE(7)=6
        NCOLTYPE(8)=6
        NCOLTYPE(9)=6
        NCOLTYPE(10)=6
        CRADIUS='m'
        LASTRADIUS=60.
        LASTWIDTH=1.0
        WINSIZE=0 !avoid compilation warning
C
        CALL BUTTON(1,'[z]oom (m)',0)
        CALL BUTTON(2,'zoom [k]',0)
        CALL BUTTON(3,'[w]hole',0)
        CALL BUTTON(4,'[s]et BG/FG',0)
        CALL BUTTON(5,'postscript',0)
        CALL BUTTON(6,'[q]uit',0)
        CALL BUTTON(7,'[i]nvert',0)
        CALL BUTTON(8,'[m]ark',0)
        CALL BUTTON(9,'[u]nmark',0)
        CALL BUTTON(10,'i[/]o marks',0)
        CALL BUTTON(11,'sa[v]e row',0)
        CALL BUTTON(12,'min[,]max',0)
        CALL BUTTON(13,'[p]anorama',0)
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
ccc10      NC1=1
        NC1=1
        NC2=NCHAN
        NS1=1
        NS2=NSCAN
        BG=A(NC1,NS1)
        FG=BG
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            BG=AMIN1(BG,A(JJ,II))
            FG=AMAX1(FG,A(JJ,II))
          END DO
        END DO
        WRITE(*,100)'Background: '
        WRITE(*,*)BG
        WRITE(*,100)'Foreground: '
        WRITE(*,*)FG
C
16      XMIN=REAL(NC1)-.6
        XMAX=REAL(NC2)+.6
        YMIN=REAL(NS1)-.6
        YMAX=REAL(NS2)+.6
        XMINA=XMIN*1.7*REAL(NBIN)
        XMAXA=XMAX*1.7*REAL(NBIN)
        YMINA=YMIN*1.7*REAL(NBIN)
        YMAXA=YMAX*1.7*REAL(NBIN)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.84)
            CALL RPGENV(XMINA,XMAXA,YMINA,YMAXA,1,-2)
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBOX('ICTMS',0.,0,'ICTMS',0.,0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,1,-2)
            CALL PGBOX('IBTNS',0.,0,'IBTNS',0.,0)
            CALL PGPTEXT(XMAX+(XMAX-XMIN)*0.10,(YMIN+YMAX)/2.,
     +       90.,0.5,'(arcsec)')
          ELSE
            CALL PGENV(XMINA,XMAXA,YMINA,YMAXA,1,-2)
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBOX('ICTMS',0.,0,'ICTMS',0.,0)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
            CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            CALL PGBOX('IBTNS',0.,0,'IBTNS',0.,0)
            CALL PGPTEXT(XMAX+(XMAX-XMIN)*0.10,(YMIN+YMAX)/2.,
     +       90.,0.5,'(arcsec)')
          END  IF
          CALL PGLABEL('pixel','pixel',FILENAME)
        END DO
        MINVAL=A(NC1,NS1)
        MAXVAL=MINVAL
        IIMAX=NS1
        IIMIN=NS1
        JJMAX=NC1
        JJMIN=NC1
        MEANVAL=0.D0
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            IF(A(JJ,II).LT.MINVAL)THEN
              MINVAL=A(JJ,II)
              IIMIN=II
              JJMIN=JJ
            END IF
            IF(A(JJ,II).GT.MAXVAL)THEN
              MAXVAL=A(JJ,II)
              IIMAX=II
              JJMAX=JJ
            END IF
            MEANVAL=MEANVAL+DBLE(A(JJ,II))
          END DO
        END DO
        NTPT=(NS2-NS1+1)*(NC2-NC1+1)
        MEANVAL=MEANVAL/DBLE(NTPT)
        WRITE(*,'(A,I5,A,I5)')'> From Scan    #',NS1,' to ',NS2
        WRITE(*,'(A,I5,A,I5)')'> From Channel #',NC1,' to ',NC2
        WRITE(*,'(A,I10)')'> Total number of pixels: ',NTPT 
        WRITE(*,100)'> Minimum: '
        WRITE(*,*)MINVAL
        WRITE(*,100)' --> in pixel:'
        WRITE(*,*)JJMIN,IIMIN
        WRITE(*,100)'> Maximum: '
        WRITE(*,*)MAXVAL
        WRITE(*,100)' --> in pixel:'
        WRITE(*,*)JJMAX,IIMAX
        WRITE(*,100)'> Mean   : '
        WRITE(*,*)REAL(MEANVAL)
        SIGMA=0.D0
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            SIGMA=SIGMA+(DBLE(A(JJ,II))-MEANVAL)*
     +       (DBLE(A(JJ,II))-MEANVAL) 
          END DO
        END DO
        SIGMA=DSQRT(SIGMA/DBLE(NTPT))
        WRITE(*,100)'> Sigma  : '
        WRITE(*,*)REAL(SIGMA)
        WRITE(*,*)
17      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGGRAY(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          IF(NOBJ.GT.0)THEN
            CALL PLOTSYMB(0,.TRUE.)
          END IF
        END DO
C------------------------------------------------------------------------------
C
18      CONTINUE
        NB=0
C------------------------------------------------------------------------------
20      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        NBLOCAL=INDEX('zkws qimu/v,p',CH)
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.1)THEN
          CALL BUTTON(1,'zoom [m]',5)
          WRITE(*,101)'Press cursor at two corners of the imaginary '
     +     //'BOX to be zoomed'
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(1,'zoom [m]',0)
            GOTO 18
          END IF
          IXC1=INT(XC+0.5)
          IYC1=INT(YC+0.5)
          IF(IXC1.LT.NC1) IXC1=NC1
          IF(IXC1.GT.NC2) IXC1=NC2
          IF(IYC1.LT.NS1) IYC1=NS1
          IF(IYC1.GT.NS2) IYC1=NS2
          WRITE(*,110)'Cursor at ',IXC1,IYC1
C
          IF(LCOLOR(1)) CALL PGSCI(5)
          CALL RPGBAND(2,0,REAL(IXC1),REAL(IYC1),XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(1,'zoom [m]',0)
            GOTO 18
          END IF
          IXC2=INT(XC+0.5)
          IYC2=INT(YC+0.5)
          IF(IXC2.LT.NC1) IXC2=NC1
          IF(IXC2.GT.NC2) IXC2=NC2
          IF(IYC2.LT.NS1) IYC2=NS1
          IF(IYC2.GT.NS2) IYC2=NS2
          WRITE(*,110)'Cursor at ',IXC2,IYC2
C
          NC1=MIN0(IXC1,IXC2)
          NC2=MAX0(IXC1,IXC2)
          NS1=MIN0(IYC1,IYC2)
          NS2=MAX0(IYC1,IYC2)
          CALL BUTTON(1,'zoom [m]',0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,'zoom [k]',5)
          WRITE(*,101)'Enter the coordinates of two corners of the '//
     +     'imaginary box to be zoomed.'
          WRITE(*,100)'First point  (channel,scan) '
          WRITE(CDUMMY,'(I4,A1,I4)')NC1,',',NS1
          CALL READ2I(CDUMMY,IXC1,IYC1)
          IF(.NOT.INSIDE1(IXC1,1,NCHAN,IYC1,1,NSCAN))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(2,'zoom [k]',0)
            GOTO 18
          END IF
          WRITE(*,100)'Second point (channel,scan) '
          WRITE(CDUMMY,'(I4,A1,I4)')NC2,',',NS2
          CALL READ2I(CDUMMY,IXC2,IYC2)
          IF(.NOT.INSIDE1(IXC2,1,NCHAN,IYC2,1,NSCAN))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(2,'zoom [k]',0)
            GOTO 18
          END IF
          NC1=MIN0(IXC1,IXC2)
          NC2=MAX0(IXC1,IXC2)
          NS1=MIN0(IYC1,IYC2)
          NS2=MAX0(IYC1,IYC2)
          CALL BUTTON(2,'zoom [k]',0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,'[w]hole',5)
          NC1=1
          NC2=NCHAN
          NS1=1
          NS2=NSCAN
          CALL BUTTON(3,'[w]hole',0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,'[s]et BG/FG',5)
          WRITE(*,100)'BACKGROUND    : '
          WRITE(*,*)BG
          WRITE(*,100)'FOREGROUND    : '
          WRITE(*,*)FG
          WRITE(*,100)'NEW BACKGROUND '
          WRITE(CDUMMY,*)BG
          BG=READF(CDUMMY)
          WRITE(*,100)'NEW FOREGROUND '
          WRITE(CDUMMY,*)FG
          FG=READF(CDUMMY)
          CALL BUTTON(4,'[s]et BG/FG',0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.5)THEN
          CALL PGQID(IDOLD)
          IDNEW=PGOPEN('?')
          WRITE(*,100) 'Creating new plot..'
          CALL PGENV(XMINA,XMAXA,YMINA,YMAXA,1,-2)
          CALL PGSCI(3)
          CALL PGBOX('ICTMS',0.,0,'ICTMS',0.,0)
          CALL PGSCI(1)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGBOX('IBTNS',0.,0,'IBTNS',0.,0)
          CALL PGPTEXT(XMAX+(XMAX-XMIN)*0.10,(YMIN+YMAX)/2.,
     +     90.,0.5,'(arcsec)')
          CALL PGLABEL('pixel','pixel',FILENAME)
          CALL PGGRAY(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          IF(NOBJ.GT.0)THEN
            CALL PLOTSYMB(0,.TRUE.)
          END IF
          CALL PGCLOS(IDNEW)
          WRITE(*,101) '   ...OK! Plot created and closed'
          CALL PGSLCT(IDOLD)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,'[q]uit',5)
          WRITE(*,100)'Do you really want to QUIT (y/n) '
          CQUIT(1:1)=READC('y','yn')
          IF(CQUIT.EQ.'y')THEN
            CALL PGEND
            STOP
          END IF
          CALL BUTTON(6,'[q]uit',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,'[i]nvert',5)
          FDUM=BG
          BG=FG
          FG=FDUM
          CALL BUTTON(7,'[i]nvert',0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,'[m]ark',5)
C
          WRITE(*,101) '[c]ircle'
          WRITE(*,101) '[s]quare'
          WRITE(*,101) '[t]riangle'
          WRITE(*,101) 'slit [1] (1 object)'
          WRITE(*,101) 'slit [2] (2 objects)'
          WRITE(*,100) 'Option (c/s/t/1/2) '
          CH(1:1)=READC('s','cst12')
          IF(CH.EQ.'s')THEN
            NTYPE(NOBJ+1)=1                                     !square
          ELSEIF(CH.EQ.'c')THEN
            NTYPE(NOBJ+1)=2                                     !circle
          ELSEIF(CH.EQ.'t')THEN
            NTYPE(NOBJ+1)=3                                     !triangle
          ELSEIF(CH.EQ.'1')THEN
            NTYPE(NOBJ+1)=4                                     !slit 1
          ELSEIF(CH.EQ.'2')THEN
            NTYPE(NOBJ+1)=5                                     !slit 2
          END IF
C
          IF(NTYPE(NOBJ+1).EQ.5)THEN
            WRITE(*,100)'Press mouse button in object #1...'
            CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
            IF(.NOT.INSIDE1(NINT(XC),NC1,NC2,NINT(YC),NS1,NS2))THEN
              WRITE(*,101)'ERROR: limits out of plot.'
              CALL BUTTON(8,'[m]ark',0)
              GOTO 18
            END IF
            XOBJ_1=XC
            YOBJ_1=YC
            WRITE(*,101)'  ..OK!'
            CALL COOR_AD(XC,YC,ARH,ARM,ARS,DECSIG,DECD,DECM,DECS)
            ARHOBJ_1=ARH
            ARMOBJ_1=ARM
            ARSOBJ_1=ARS
            DECSIGOBJ_1=DECSIG
            DECDOBJ_1=DECD
            DECMOBJ_1=DECM
            DECSOBJ_1=DECS
            WRITE(*,100)'Press mouse button in object #2...'
            CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
            IF(.NOT.INSIDE1(NINT(XC),NC1,NC2,NINT(YC),NS1,NS2))THEN
              WRITE(*,101)'ERROR: limits out of plot.'
              CALL BUTTON(8,'[m]ark',0)
              GOTO 18
            END IF
            XOBJ_2=XC
            YOBJ_2=YC
            WRITE(*,101)'  ..OK!'
            CALL COOR_AD(XC,YC,ARH,ARM,ARS,DECSIG,DECD,DECM,DECS)
            ARHOBJ_2=ARH
            ARMOBJ_2=ARM
            ARSOBJ_2=ARS
            DECSIGOBJ_2=DECSIG
            DECDOBJ_2=DECD
            DECMOBJ_2=DECM
            DECSOBJ_2=DECS
            CALL PACENT(ARHOBJ_1,ARMOBJ_1,ARSOBJ_1,
     +                  DECSIGOBJ_1,DECDOBJ_1,DECMOBJ_1,DECSOBJ_1,
     +                  ARHOBJ_2,ARMOBJ_2,ARSOBJ_2,
     +                  DECSIGOBJ_2,DECDOBJ_2,DECMOBJ_2,DECSOBJ_2,
     +                  ARH,ARM,ARS,
     +                  DECSIG,DECD,DECM,DECS,
     +                  PA_12)
            ARHOBJ(NOBJ+1)=ARH
            ARMOBJ(NOBJ+1)=ARM
            ARSOBJ(NOBJ+1)=ARS
            DECSIGOBJ(NOBJ+1)=DECSIG
            DECDOBJ(NOBJ+1)=DECD
            DECMOBJ(NOBJ+1)=DECM
            DECSOBJ(NOBJ+1)=DECS
            PAOBJ(NOBJ+1)=PA_12
            XOBJ(NOBJ+1)=(XOBJ_1+XOBJ_2)/2.
            YOBJ(NOBJ+1)=(YOBJ_1+YOBJ_2)/2.
            WRITE(*,101) 'Centroid:'
            CALL COOR_AD(XOBJ(NOBJ+1),YOBJ(NOBJ+1),
     +       ARH,ARM,ARS,DECSIG,DECD,DECM,DECS)
            WRITE(*,100) 'Slit P.A. (deg.): '
            WRITE(*,*) PAOBJ(NOBJ+1)
          ELSE
            WRITE(*,100)'Press mouse button in object...'
            CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
            IF(.NOT.INSIDE1(NINT(XC),NC1,NC2,NINT(YC),NS1,NS2))THEN
              WRITE(*,101)'ERROR: limits out of plot.'
              CALL BUTTON(8,'[m]ark',0)
              GOTO 18
            END IF
            XOBJ(NOBJ+1)=XC
            YOBJ(NOBJ+1)=YC
            WRITE(*,101)'  ..OK!'
            CALL COOR_AD(XC,YC,ARH,ARM,ARS,DECSIG,DECD,DECM,DECS)
            ARHOBJ(NOBJ+1)=ARH
            ARMOBJ(NOBJ+1)=ARM
            ARSOBJ(NOBJ+1)=ARS
            DECSIGOBJ(NOBJ+1)=DECSIG
            DECDOBJ(NOBJ+1)=DECD
            DECMOBJ(NOBJ+1)=DECM
            DECSOBJ(NOBJ+1)=DECS
          END IF
C
          IF((NTYPE(NOBJ+1).GE.1).AND.(NTYPE(NOBJ+1).LE.3))THEN
            WRITE(*,100)'Radius with [m]ouse or [k]eyboard (m/s) '
            CRADIUS(1:1)=READC(CRADIUS,'mk')
            IF(CRADIUS.EQ.'m')THEN
              WRITE(*,100)'Press mouse button in radius...'
              IF(LCOLOR(1)) CALL PGSCI(5)
              CALL RPGBAND(1,0,XC,YC,XC,YC,CH)
              IF(LCOLOR(1)) CALL PGSCI(1)
              IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
                WRITE(*,101)'ERROR: mouse buttom has not been detected.'
                CALL BUTTON(8,'[m]ark',0)
                GOTO 18
              END IF
              IF(.NOT.INSIDE1(NINT(XC),NC1,NC2,NINT(YC),NS1,NS2))THEN
                WRITE(*,101)'ERROR: limits out of plot.'
                CALL BUTTON(8,'[m]ark',0)
                GOTO 18
              END IF
              WRITE(*,101)'  ..OK!'
              ROBJ(NOBJ+1)=SQRT((XOBJ(NOBJ+1)-XC)*(XOBJ(NOBJ+1)-XC)+
     +         (YOBJ(NOBJ+1)-YC)*(YOBJ(NOBJ+1)-YC))
              RARCOBJ(NOBJ+1)=ROBJ(NOBJ+1)*1.7*REAL(NBIN)
            ELSE
              WRITE(CDUMMY,*)LASTRADIUS
              WRITE(*,100)'Radius (arcsec) '
              RARCOBJ(NOBJ+1)=READF(CDUMMY)
              LASTRADIUS=RARCOBJ(NOBJ+1)
              ROBJ(NOBJ+1)=RARCOBJ(NOBJ+1)/(1.7*REAL(NBIN))
            END IF
            WRITE(*,100)'Text description (max. 20 char.) '
            READ(*,'(A)')TEXTOBJ(NOBJ+1)
          ELSEIF(NTYPE(NOBJ+1).EQ.4)THEN
            WRITE(CDUMMY,*)LASTRADIUS
            WRITE(*,100)'Slit length (arcsec) '
            RARCOBJ(NOBJ+1)=READF(CDUMMY)
            LASTRADIUS=RARCOBJ(NOBJ+1)
            ROBJ(NOBJ+1)=RARCOBJ(NOBJ+1)/(1.7*REAL(NBIN))
            WRITE(*,100)'Position Angle (degrees)'
            PAOBJ(NOBJ+1)=READF('@')
            WRITE(*,100)'Slit width (arcsec) '
            SWARCOBJ(NOBJ+1)=READF('@')
            SWOBJ(NOBJ+1)=SWARCOBJ(NOBJ+1)/(1.7*REAL(NBIN))
            WRITE(TEXTOBJ(NOBJ+1),137) PAOBJ(NOBJ+1),SWARCOBJ(NOBJ+1)
          ELSEIF(NTYPE(NOBJ+1).EQ.5)THEN
            WRITE(CDUMMY,*)LASTRADIUS
            WRITE(*,100)'Slit length (arcsec) '
            RARCOBJ(NOBJ+1)=READF(CDUMMY)
            LASTRADIUS=RARCOBJ(NOBJ+1)
            ROBJ(NOBJ+1)=RARCOBJ(NOBJ+1)/(1.7*REAL(NBIN))
            WRITE(CDUMMY,*)LASTWIDTH
            WRITE(*,100)'Slit width (arcsec) '
            SWARCOBJ(NOBJ+1)=READF(CDUMMY)
            LASTWIDTH=SWARCOBJ(NOBJ+1)
            SWOBJ(NOBJ+1)=SWARCOBJ(NOBJ+1)/(1.7*REAL(NBIN))
            WRITE(TEXTOBJ(NOBJ+1),137) PAOBJ(NOBJ+1),SWARCOBJ(NOBJ+1)
          END IF
C
          NOBJ=NOBJ+1
C
          N0=NOBJ
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PLOTSYMB(N0,.FALSE.)
          END DO
C
          CALL BUTTON(8,'[m]ark',0)
          GOTO 18
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.9)THEN
          CALL BUTTON(9,'[u]nmark',5)
          IF(NOBJ.EQ.0) GOTO 18
          IF(NOBJ.EQ.1)THEN
            NOBJ=0
            CALL BUTTON(9,'[u]nmark',0)
            GOTO 17
          END IF
C
          WRITE(*,100)'Press mouse button...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(9,'[u]nmark',0)
            GOTO 18
          END IF
          IF(.NOT.INSIDE1(NINT(XC),NC1,NC2,NINT(YC),NS1,NS2))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(9,'[u]nmark',0)
            GOTO 18
          END IF
          WRITE(*,101)'  ..OK!'
          CALL REMOVESYMB(XC,YC,.TRUE.)
C
          CALL BUTTON(9,'[u]nmark',0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.10)THEN
          CALL BUTTON(10,'i[/]o marks',5)
          WRITE(*,100)'[l]ist, [i]nput or [o]utput mark file, or '//
     +     ' [c]olors (l,i/o/c) '
          CIOMARKS(1:1)=READC('@','lioc')
          IF(CIOMARKS.EQ.'l')THEN
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,101)'HH MM SS.dd +DD MM SS.d RRRR.dd TT TEXT(*20)'
            DO N0=1,NOBJ
              WRITE(*,333)ARHOBJ(N0),ARMOBJ(N0),ARSOBJ(N0),
     +         DECSIGOBJ(N0),DECDOBJ(N0),DECMOBJ(N0),DECSOBJ(N0),
     +         RARCOBJ(N0),NTYPE(N0),TEXTOBJ(N0)
            END DO
            WRITE(*,101)'--------------------------------------------'
            WRITE(*,*)
          ELSEIF(CIOMARKS.EQ.'i')THEN
            WRITE(*,101)'--------------------------------------------'//
     +       '-----------'
            WRITE(*,101)'Required file format: '
            WRITE(*,101)'--------------------------------------------'//
     +       '-----------'
            WRITE(*,101)'00000000011111111112222222222333333333344444'//
     +       '44444555555'
            WRITE(*,101)'12345678901234567890123456789012345678901234'//
     +       '56789012345'
            WRITE(*,101)'HH MM SS.dd +DD MM SS.d RRRR.dd TT TEXT(*20)'
            WRITE(*,*)
            WRITE(*,101)' RRRR.dd - symbol radius (arcsecs)'
            WRITE(*,101)' TT      - symbol type: 1-square'
            WRITE(*,101)'                        2-circle'
            WRITE(*,101)'                        3-triangle'
            WRITE(*,101)'                        4-slit #1'
            WRITE(*,101)'                        5-slit #2'
            WRITE(*,101)'--------------------------------------------'//
     +       '-----------'
            WRITE(*,100)'Marks file name'
            MARKFILE(1:75)=READC('@','@')
            OPEN(77,FILE=MARKFILE,STATUS='OLD',FORM='FORMATTED')
            IF(NOBJ.GT.0)THEN
              WRITE(*,101)' [a]ppend new marks to existing marks'
              WRITE(*,101)' [d]elete previous marks and use new ones'
              WRITE(*,100)'Option (a/d) '
              CMERGE(1:1)=READC('@','ad')
              IF(CMERGE.EQ.'a')THEN
              ELSEIF(CMERGE.EQ.'d')THEN
                NOBJ=0
              END IF
            END IF
            WRITE(*,100)'J2000.0 (y/n) '
            CJ2000(1:1)=READC('y','yn')
            IF(CJ2000.EQ.'n')THEN
              WRITE(*,100)'Equinox'
              EQUINOXFILE=READF('@')
            END IF
C
            READ(77,*)                              !ignoramos la primera linea
300         READ(77,333,END=310)ARHOBJ(NOBJ+1),ARMOBJ(NOBJ+1),
     +       ARSOBJ(NOBJ+1),DECSIGOBJ(NOBJ+1),DECDOBJ(NOBJ+1),
     +       DECMOBJ(NOBJ+1),DECSOBJ(NOBJ+1),RARCOBJ(NOBJ+1),
     +       NTYPE(NOBJ+1),TEXTOBJ(NOBJ+1)
            NOBJ=NOBJ+1
            WRITE(*,333)ARHOBJ(NOBJ),ARMOBJ(NOBJ),
     +       ARSOBJ(NOBJ),DECSIGOBJ(NOBJ),DECDOBJ(NOBJ),
     +       DECMOBJ(NOBJ),DECSOBJ(NOBJ),RARCOBJ(NOBJ),
     +       NTYPE(NOBJ),TEXTOBJ(NOBJ)
            IF(CJ2000.EQ.'n')THEN
              CALL PRECESOBJ(EQUINOXFILE,ARHOBJ(NOBJ),ARMOBJ(NOBJ),
     +         ARSOBJ(NOBJ),DECSIGOBJ(NOBJ),DECDOBJ(NOBJ),
     +         DECMOBJ(NOBJ),DECSOBJ(NOBJ))
              WRITE(*,334)ARHOBJ(NOBJ),ARMOBJ(NOBJ),
     +         ARSOBJ(NOBJ),DECSIGOBJ(NOBJ),DECDOBJ(NOBJ),
     +         DECMOBJ(NOBJ),DECSOBJ(NOBJ),'(J2000.00)'
            END IF
            ROBJ(NOBJ)=RARCOBJ(NOBJ)/(1.7*REAL(NBIN))          !radio en pixels
            IF((NTYPE(NOBJ).EQ.4).OR.(NTYPE(NOBJ).EQ.5))THEN    !es una rendija
              READ(TEXTOBJ(NOBJ),137)PAOBJ(NOBJ),SWARCOBJ(NOBJ)
              SWOBJ(NOBJ)=SWARCOBJ(NOBJ)/(1.7*REAL(NBIN))
            END IF
            IF(DECSIGOBJ(NOBJ+1).EQ.' ') DECSIGOBJ(NOBJ+1)='+'
            CALL COOR_XY(ARHOBJ(NOBJ),ARMOBJ(NOBJ),ARSOBJ(NOBJ),
     +       DECSIGOBJ(NOBJ),DECDOBJ(NOBJ),DECMOBJ(NOBJ),DECSOBJ(NOBJ),
     +       XOBJ(NOBJ),YOBJ(NOBJ))
            GOTO 300
310         CLOSE(77)
          ELSEIF(CIOMARKS.EQ.'o')THEN
            WRITE(*,100)'Marks file name? '
            MARKFILE(1:75)=READC('@','@')
            INQUIRE(FILE=MARKFILE,EXIST=LOGFILE)
            IF(LOGFILE)THEN
              WRITE(*,101)'ERROR: this file already exist.'
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
            ELSE
              OPEN(77,FILE=MARKFILE,STATUS='NEW',FORM='FORMATTED')
              WRITE(*,*)
              WRITE(*,101)'*** WARNING ***'
              WRITE(*,101)'This file will save J2000.0 coordinates.'
              WRITE(*,*)
              WRITE(*,100)'(press <RETURN> to continue)'
              READ(*,*)
              WRITE(*,*)
              WRITE(77,101)'HH MM SS.dd +DD MM SS.d RRRR.dd'//
     +         ' TT TEXT(*20)'
              DO N0=1,NOBJ
                WRITE(77,333)ARHOBJ(N0),ARMOBJ(N0),ARSOBJ(N0),
     +           DECSIGOBJ(N0),DECDOBJ(N0),DECMOBJ(N0),DECSOBJ(N0),
     +           RARCOBJ(N0),NTYPE(N0),TEXTOBJ(N0)
              END DO
              CLOSE(77)
            END IF
          ELSEIF(CIOMARKS.EQ.'c')THEN
            WRITE(*,100)'Color for marks type  1 '
            WRITE(CDUMMY,*) NCOLTYPE(1)
            NCOLTYPE(1)=READI(CDUMMY)
            WRITE(*,100)'Color for marks type  2 '
            WRITE(CDUMMY,*) NCOLTYPE(2)
            NCOLTYPE(2)=READI(CDUMMY)
            WRITE(*,100)'Color for marks type  3 '
            WRITE(CDUMMY,*) NCOLTYPE(3)
            NCOLTYPE(3)=READI(CDUMMY)
            WRITE(*,100)'Color for marks type  4 '
            WRITE(CDUMMY,*) NCOLTYPE(4)
            NCOLTYPE(4)=READI(CDUMMY)
          END IF
          CALL BUTTON(10,'i[/]o marks',0)
          IF(CIOMARKS.EQ.'i')THEN
            GOTO 17
          ELSE
            GOTO 18
          END IF
333       FORMAT(2(I2.2,1X),F5.2,1X,A1,2(I2.2,1X),F4.1,1X,F7.2,1X,
     +     I2.2,1X,A20)
334       FORMAT(2(I2.2,1X),F5.2,1X,A1,2(I2.2,1X),F4.1,1X,A10)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(11,'sa[v]e row',5)
          WRITE(*,101)'NOTE: displayed region will be saved.'
          WRITE(*,100)'Output file name? '
          OUTFILE(1:75)=READC('@','@')
          OPEN(34,FILE=OUTFILE,STATUS='NEW',FORM='UNFORMATTED')
          WRITE(*,100)'Saving file...'
          DO II=NS1,NS2
            WRITE(34) (A(JJ,II),JJ=NC1,NC2)
          END DO
          CLOSE(34)
          WRITE(*,101)'  ..OK!'
          WRITE(*,110)'No. of scans saved...: ',NS2-NS1+1
          WRITE(*,110)'No. of channels saved: ',NC2-NC1+1
          CALL BUTTON(11,'sa[v]e row',0)
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(12,'min[,]max',5)
          BG=REAL(MINVAL)
          FG=REAL(MAXVAL)
          CALL BUTTON(12,'min[,]max',0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.13)THEN
          CALL BUTTON(13,'[p]anorama',5)
          IF(NOBJ.EQ.0)THEN
            WRITE(*,101)'WARNING: No. of marks = 0!'
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
          ELSE
            WRITE(*,100)'Window size: [a]utomatic or [f]ixed (a/f) '
            CWSIZE(1:1)=READC('a','af')
            IF(CWSIZE.EQ.'f')THEN
              WRITE(*,100)'Window side (arcsec) '
              WINSIZE=READF('@')
              WINSIZE=WINSIZE/(1.7*REAL(NBIN))                       !en pixels
            END IF
            WRITE(*,100)'Mark type 1,2,3 or [a]ll (1/2/3/a) '
            CSELMARK(1:1)=READC('a','123a')
            DO N0=1,NOBJ
              IF(
     +         (CSELMARK.EQ.'a').OR.
     +         ((CSELMARK.EQ.'1').AND.(NTYPE(N0).EQ.1)).OR.
     +         ((CSELMARK.EQ.'2').AND.(NTYPE(N0).EQ.2)).OR.
     +         ((CSELMARK.EQ.'3').AND.(NTYPE(N0).EQ.3))
     +         )THEN
                IF(CWSIZE.EQ.'a') WINSIZE=4.*ROBJ(N0)
                FNC1=XOBJ(N0)-WINSIZE/2.
                IF(FNC1.LT.1.) FNC1=1.
                FNC2=XOBJ(N0)+WINSIZE/2.
                IF(FNC2.GT.REAL(NCHAN)) FNC2=REAL(NCHAN)
                FNS1=YOBJ(N0)-WINSIZE/2.
                IF(FNS1.LT.1.) FNS1=1.
                FNS2=YOBJ(N0)+WINSIZE/2.
                IF(FNS2.GT.REAL(NSCAN)) FNS2=REAL(NSCAN)
                XMIN=FNC1-.6
                XMAX=FNC2+.6
                YMIN=FNS1-.6
                YMAX=FNS2+.6
                XMINA=XMIN*1.7*REAL(NBIN)
                XMAXA=XMAX*1.7*REAL(NBIN)
                YMINA=YMIN*1.7*REAL(NBIN)
                YMAXA=YMAX*1.7*REAL(NBIN)
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(ITERM.EQ.1)THEN
                    CALL RPGERASW(0.,1.,0.,0.84)
                    CALL RPGENV(XMINA,XMAXA,YMINA,YMAXA,1,-2)
                    IF(LCOLOR(ITERM)) CALL PGSCI(3)
                    CALL PGBOX('ICTMS',0.,0,'ICTMS',0.,0)
                    IF(LCOLOR(ITERM)) CALL PGSCI(1)
                    CALL RPGENV(XMIN,XMAX,YMIN,YMAX,1,-2)
                    CALL PGBOX('IBTNS',0.,0,'IBTNS',0.,0)
                    CALL PGPTEXT(XMAX+(XMAX-XMIN)*0.10,(YMIN+YMAX)/2.,
     +               90.,0.5,'(arcsec)')
                  ELSE
                    CALL PGENV(XMINA,XMAXA,YMINA,YMAXA,1,-2)
                    IF(LCOLOR(ITERM)) CALL PGSCI(3)
                    CALL PGBOX('ICTMS',0.,0,'ICTMS',0.,0)
                    IF(LCOLOR(ITERM)) CALL PGSCI(1)
                    CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
                    CALL PGBOX('IBTNS',0.,0,'IBTNS',0.,0)
                    CALL PGPTEXT(XMAX+(XMAX-XMIN)*0.10,(YMIN+YMAX)/2.,
     +               90.,0.5,'(arcsec)')
                  END IF
                  CALL PGLABEL('channel','scan',' ')
                END DO
                NC1=NINT(FNC1)-.6
                IF(NC1.LT.1) NC1=1
                NC2=NINT(FNC2)+.6
                IF(NC2.GT.NCHAN) NC2=NCHAN
                NS1=NINT(FNS1)-.6
                IF(NS1.LT.1) NS1=1
                NS2=NINT(FNS2)+.6
                IF(NS2.GT.NSCAN) NS2=NSCAN
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGGRAY(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
                  CALL PLOTSYMB(N0,.FALSE.)
                END DO
                WRITE(*,*)
                WRITE(*,101)'>>> Object: '//TEXTOBJ(N0)
                WRITE(*,100)'Next object (y/n) '
                CNEXT(1:1)=READC('y','yn')
                IF(CNEXT.EQ.'n')THEN
                  NC1=1
                  NC2=NCHAN
                  NS1=1
                  NS2=NSCAN
                  CALL BUTTON(13,'[p]anorama',0)
                  GOTO 16
                END IF
              END IF
            END DO
            NC1=1
            NC2=NCHAN
            NS1=1
            NS2=NSCAN
            CALL BUTTON(13,'[p]anorama',0)
            GOTO 16
          END IF
          CALL BUTTON(13,'[p]anorama',0)
C------------------------------------------------------------------------------
        ELSE
          JJ=INT(XC+0.5)
          II=INT(YC+0.5)
          IF(INSIDE1(JJ,NC1,NC2,II,NS1,NS2))THEN
            WRITE(*,111)'Cursor at ',XC,YC,'      Pixel value: '
            WRITE(*,*)A(JJ,II)
            CALL COOR_AD(XC,YC,ARH,ARM,ARS,DECSIG,DECD,DECM,DECS)
            IF(NOBJ.GT.0) CALL REMOVESYMB(XC,YC,.FALSE.)
            WRITE(*,*)
          ELSE
            WRITE(*,101)'ERROR: Cursor out of plot.'
          END IF
        END IF
C------------------------------------------------------------------------------
        GOTO 20
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I5,2X,I5)
111     FORMAT(A,F8.2,2X,F8.2,A,$)
137     FORMAT(F7.2,2X,F6.2)
C
        END
C
C**********************************************************************
C
C Si el punto J0,I0 esta dentro del rectangulo definido por J1,J2,I1,I2
C la funcion devuelve .TRUE.
        LOGICAL FUNCTION INSIDE1(J0,J1,J2,I0,I1,I2)
        IMPLICIT NONE
        INTEGER J0,J1,J2
        INTEGER I0,I1,I2
C
        INSIDE1=.TRUE.
        IF(J0.LT.J1) GOTO 70
        IF(J0.GT.J2) GOTO 70
        IF(I0.LT.I1) GOTO 70
        IF(I0.GT.I2) GOTO 70
        RETURN
70      CONTINUE
        INSIDE1=.FALSE.
        RETURN
        END
C
C******************************************************************************
C
C Dibuja un simbolo centrado en XOBJ,YOBJ, con radio ROBJ en los NOBJ objectos
C introducidos. Si N0=0 dibujamos todos los simbolos definidos. En caso
C contrario, solo pintamos el simbolo N0. Para las rendijas dibujamos con el
C angulo de posicion dado por PAOBJ, y una anchura dada por SWOBJ
C LINFORM=.TRUE. informa de los objetos que caen fuera
        SUBROUTINE PLOTSYMB(N0,LINFORM)
        IMPLICIT NONE
        INTEGER TRUELEN
C
        INTEGER N0
C
        INTEGER NMAXOBJ
        PARAMETER (NMAXOBJ=8000)
        REAL PI
        PARAMETER(PI=3.141592654)
C
        INTEGER NOBJ
        INTEGER NTYPE(NMAXOBJ)
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NCOLTYPE(10)
ccc        INTEGER NTERM,IDN(8)
        REAL XOBJ(NMAXOBJ),YOBJ(NMAXOBJ)
        REAL ROBJ(NMAXOBJ),PAOBJ(NMAXOBJ),SWOBJ(NMAXOBJ)
        CHARACTER*20 TEXTOBJ(NMAXOBJ)
ccc        LOGICAL LCOLOR(8)
C
        COMMON/BLKSYMB1/NOBJ
        COMMON/BLKSYMB2/XOBJ,YOBJ,ROBJ,PAOBJ,SWOBJ
        COMMON/BLKSYMB3/NTYPE
        COMMON/BLKSYMB7/TEXTOBJ
        COMMON/BLKSYMB8/NC1,NC2,NS1,NS2
        COMMON/BLKSYMB9/NCOLTYPE
ccc        COMMON/BLKDEVICE1/NTERM,IDN
ccc        COMMON/BLKDEVICE2/LCOLOR
C
        INTEGER N1,N2,N,I
        INTEGER LW
        REAL X(37),Y(37),ANG
        LOGICAL LINFORM
        LOGICAL INSIDE1
C------------------------------------------------------------------------------
        IF((N0.LT.0).OR.(N0.GT.NOBJ))THEN
          WRITE(*,101)'ERROR: symbol number out of range.'
          RETURN
        END IF
        IF(N0.EQ.0)THEN
          N1=1
          N2=NOBJ
        ELSE
          N1=N0
          N2=N0
        END IF
C
        DO N=N1,N2
          CALL PGSCI(NCOLTYPE(NTYPE(N)))
          IF(NTYPE(N).EQ.2)THEN                                       !circulos
            DO I=1,37
              ANG=REAL(I-1)*3.141593/18.
              X(I)=ROBJ(N)*COS(ANG)+XOBJ(N)
              Y(I)=ROBJ(N)*SIN(ANG)+YOBJ(N)
            END DO
            CALL PGLINE(37,X,Y)
          ELSEIF(NTYPE(N).EQ.1)THEN                                  !cuadrados
            CALL PGSFS(2)
            CALL PGRECT(XOBJ(N)-ROBJ(N),XOBJ(N)+ROBJ(N),
     +       YOBJ(N)-ROBJ(N),YOBJ(N)+ROBJ(N))
            CALL PGSFS(1)
          ELSEIF(NTYPE(N).EQ.3)THEN                                 !triangulos
            X(1)=XOBJ(N)
            Y(1)=YOBJ(N)+ROBJ(N)
            X(2)=XOBJ(N)-ROBJ(N)*.8660
            Y(2)=YOBJ(N)-ROBJ(N)*.5000
            X(3)=XOBJ(N)+ROBJ(N)*.8660
            Y(3)=YOBJ(N)-ROBJ(N)*.5000
            X(4)=XOBJ(N)
            Y(4)=YOBJ(N)+ROBJ(N)
            CALL PGLINE(4,X,Y)
          ELSEIF(NTYPE(N).EQ.4)THEN                             !rendija tipo 1
            X(1)=XOBJ(N)+0.5*SWOBJ(N)*COS(PAOBJ(N)*PI/180.)
            Y(1)=YOBJ(N)+0.5*SWOBJ(N)*SIN(PAOBJ(N)*PI/180.)
            X(2)=X(1)+ROBJ(N)*COS((PAOBJ(N)+90.)*PI/180.)
            Y(2)=Y(1)+ROBJ(N)*SIN((PAOBJ(N)+90.)*PI/180.)
            CALL PGLINE(2,X,Y)
            X(1)=XOBJ(N)-0.5*SWOBJ(N)*COS(PAOBJ(N)*PI/180.)
            Y(1)=YOBJ(N)-0.5*SWOBJ(N)*SIN(PAOBJ(N)*PI/180.)
            X(2)=X(1)+ROBJ(N)*COS((PAOBJ(N)+90.)*PI/180.)
            Y(2)=Y(1)+ROBJ(N)*SIN((PAOBJ(N)+90.)*PI/180.)
            CALL PGLINE(2,X,Y)
          ELSEIF(NTYPE(N).EQ.5)THEN                             !rendija tipo 2
            X(1)=XOBJ(N)+0.5*SWOBJ(N)*COS(PAOBJ(N)*PI/180.)
            Y(1)=YOBJ(N)+0.5*SWOBJ(N)*SIN(PAOBJ(N)*PI/180.)
            X(2)=X(1)+ROBJ(N)/2.0*COS((PAOBJ(N)+90.)*PI/180.)
            Y(2)=Y(1)+ROBJ(N)/2.0*SIN((PAOBJ(N)+90.)*PI/180.)
            CALL PGLINE(2,X,Y)
            X(1)=XOBJ(N)-0.5*SWOBJ(N)*COS(PAOBJ(N)*PI/180.)
            Y(1)=YOBJ(N)-0.5*SWOBJ(N)*SIN(PAOBJ(N)*PI/180.)
            X(2)=X(1)+ROBJ(N)/2.0*COS((PAOBJ(N)+90.)*PI/180.)
            Y(2)=Y(1)+ROBJ(N)/2.0*SIN((PAOBJ(N)+90.)*PI/180.)
            CALL PGLINE(2,X,Y)
            X(1)=XOBJ(N)+0.5*SWOBJ(N)*COS(PAOBJ(N)*PI/180.)
            Y(1)=YOBJ(N)+0.5*SWOBJ(N)*SIN(PAOBJ(N)*PI/180.)
            X(2)=X(1)+ROBJ(N)/2.0*COS((PAOBJ(N)-90.)*PI/180.)
            Y(2)=Y(1)+ROBJ(N)/2.0*SIN((PAOBJ(N)-90.)*PI/180.)
            CALL PGLINE(2,X,Y)
            X(1)=XOBJ(N)-0.5*SWOBJ(N)*COS(PAOBJ(N)*PI/180.)
            Y(1)=YOBJ(N)-0.5*SWOBJ(N)*SIN(PAOBJ(N)*PI/180.)
            X(2)=X(1)+ROBJ(N)/2.0*COS((PAOBJ(N)-90.)*PI/180.)
            Y(2)=Y(1)+ROBJ(N)/2.0*SIN((PAOBJ(N)-90.)*PI/180.)
            CALL PGLINE(2,X,Y)
          END IF
          IF((NTYPE(N).NE.4).AND.(NTYPE(N).NE.5))THEN
            IF(INSIDE1(NINT(XOBJ(N)),NC1,NC2,NINT(YOBJ(N)),NS1,NS2))THEN
              LW=TRUELEN(TEXTOBJ(N))
              IF(LW.GT.0)THEN
                CALL PGSCH(0.6)
                CALL PGPTEXT(XOBJ(N),YOBJ(N)+ROBJ(N)*1.2,0.0,0.5,
     +           TEXTOBJ(N)(1:LW))
                CALL PGSCH(1.0)
              END IF
            ELSE
              IF(LINFORM)THEN
                WRITE(*,101)'>>> WARNING: object out of field: '//
     +           TEXTOBJ(N)
              END IF
            END IF
          END IF
          CALL PGSCI(1)
        END DO
        CALL PGSCI(1)
C
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Elimina el simbolo mas cercano a XC,YC si LREMOVE=.TRUE.
C Si LREMOVE=.FALSE. la subrutina simplemente informa acerca del objeto mas
C cercano.
        SUBROUTINE REMOVESYMB(XC,YC,LREMOVE)
        IMPLICIT NONE
        REAL XC,YC
        LOGICAL LREMOVE
C
        INTEGER NMAXOBJ
        PARAMETER (NMAXOBJ=8000)
C
        INTEGER NOBJ
        REAL XOBJ(NMAXOBJ),YOBJ(NMAXOBJ)
        REAL ROBJ(NMAXOBJ),PAOBJ(NMAXOBJ)
        REAL SWARCOBJ(NMAXOBJ),SWOBJ(NMAXOBJ)
        INTEGER NTYPE(NMAXOBJ)
        INTEGER ARHOBJ(NMAXOBJ)
        INTEGER ARMOBJ(NMAXOBJ)
        REAL ARSOBJ(NMAXOBJ)
        CHARACTER*1 DECSIGOBJ(NMAXOBJ)
        INTEGER DECDOBJ(NMAXOBJ)
        INTEGER DECMOBJ(NMAXOBJ)
        REAL DECSOBJ(NMAXOBJ)
        REAL RARCOBJ(NMAXOBJ)          !tamanho de la marca en segundos de arco
        CHARACTER*20 TEXTOBJ(NMAXOBJ)
C
        COMMON/BLKSYMB1/NOBJ
        COMMON/BLKSYMB2/XOBJ,YOBJ,ROBJ,PAOBJ,SWOBJ
        COMMON/BLKSYMB3/NTYPE
        COMMON/BLKSYMB4/ARHOBJ,ARMOBJ,ARSOBJ
        COMMON/BLKSYMB5/DECSIGOBJ,DECDOBJ,DECMOBJ,DECSOBJ
        COMMON/BLKSYMB6/RARCOBJ,SWARCOBJ
        COMMON/BLKSYMB7/TEXTOBJ
C
        INTEGER N,N0
        REAL RMIN,DIST
C------------------------------------------------------------------------------
        N0=0 !avoid compilation warning
        RMIN=1.0E10
        DO N=1,NOBJ
          DIST=(XC-XOBJ(N))*(XC-XOBJ(N))+(YC-YOBJ(N))*(YC-YOBJ(N))
          IF(DIST.LT.RMIN)THEN
            RMIN=DIST
            N0=N
          END IF
        END DO
C
        IF(.NOT.LREMOVE)THEN
          WRITE(*,101)'>>> Nearest symbol is object: '//TEXTOBJ(N0)
          RETURN
        ELSE
          WRITE(*,101)'>>> Symbol removed is object: '//TEXTOBJ(N0)
          WRITE(*,*)
        ENDIF
C N0 corresponde al simbolo mas cercano. Lo eliminamos de la lista
        IF(N0.LT.NOBJ)THEN
          DO N=N0,NOBJ-1
            XOBJ(N)=XOBJ(N+1)
            YOBJ(N)=YOBJ(N+1)
            ROBJ(N)=ROBJ(N+1)
            NTYPE(N)=NTYPE(N+1)
            ARHOBJ(N)=ARHOBJ(N+1)
            ARMOBJ(N)=ARMOBJ(N+1)
            ARSOBJ(N)=ARSOBJ(N+1)
            DECSIGOBJ(N)=DECSIGOBJ(N+1)
            DECDOBJ(N)=DECDOBJ(N+1)
            DECMOBJ(N)=DECMOBJ(N+1)
            DECSOBJ(N)=DECSOBJ(N+1)
            RARCOBJ(N)=RARCOBJ(N+1)
            PAOBJ(N)=PAOBJ(N+1)
            SWOBJ(N)=SWOBJ(N+1)
            SWARCOBJ(N)=SWARCOBJ(N+1)
            TEXTOBJ(N)=TEXTOBJ(N+1)
          END DO
        END IF
        NOBJ=NOBJ-1
C
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C
C Introducimos posicion del cursor y calcula R.A. and DEC.
C
        SUBROUTINE COOR_AD(XC,YC,ARH,ARM,ARS,DECSIG,DECD,DECM,DECS)
        IMPLICIT NONE
        REAL XC,YC
        INTEGER ARH,ARM
        REAL ARS
        CHARACTER*1 DECSIG
        INTEGER DECD,DECM
        REAL DECS
C
        DOUBLE PRECISION PI
        PARAMETER(PI=3.14159265358979323846D0)
C
        INTEGER NBIN
        INTEGER CNPIX1,CNPIX2
        INTEGER PLTRAH,PLTRAM
        INTEGER PLTDECD,PLTDECM
        DOUBLE PRECISION ARHD,ARMD,ARSD
        DOUBLE PRECISION DECDD,DECMD,DECSD
        DOUBLE PRECISION PPO3,PPO6
        DOUBLE PRECISION XPIXELSZ,YPIXELSZ
        DOUBLE PRECISION AMDX1,AMDX2,AMDX3,AMDX4,AMDX5,AMDX6,AMDX7
        DOUBLE PRECISION AMDX8,AMDX9,AMDX10,AMDX11,AMDX12,AMDX13
        DOUBLE PRECISION AMDY1,AMDY2,AMDY3,AMDY4,AMDY5,AMDY6,AMDY7
        DOUBLE PRECISION AMDY8,AMDY9,AMDY10,AMDY11,AMDY12,AMDY13
        DOUBLE PRECISION X,Y,XX,YY
        DOUBLE PRECISION RA0,DEC0,RA1,DEC1
        DOUBLE PRECISION PLTRAS,PLTDECS
        CHARACTER*1 DECSIGD
        CHARACTER*1 PLTDECSN
C
        COMMON/BLKNBIN/NBIN
        COMMON/BLKPOSS1/PPO3,PPO6,XPIXELSZ,YPIXELSZ
        COMMON/BLKPOSS2/AMDX1,AMDX2,AMDX3,AMDX4,AMDX5,AMDX6,AMDX7,
     +   AMDX8,AMDX9,AMDX10,AMDX11,AMDX12,AMDX13
        COMMON/BLKPOSS3/AMDY1,AMDY2,AMDY3,AMDY4,AMDY5,AMDY6,AMDY7,
     +   AMDY8,AMDY9,AMDY10,AMDY11,AMDY12,AMDY13
        COMMON/BLKPOSS4/CNPIX1,CNPIX2
        COMMON/BLKPOSS5/PLTRAH,PLTRAM,PLTDECD,PLTDECM
        COMMON/BLKPOSS6/PLTRAS,PLTDECS
        COMMON/BLKPOSS7/PLTDECSN
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        RA0=DBLE(PLTRAH)+DBLE(PLTRAM)/60.D0+PLTRAS/3600.D0
        RA0=RA0*15.D0
        RA0=RA0*PI/180.D0
        DEC0=DBLE(PLTDECD)+DBLE(PLTDECM)/60.D0+PLTDECS/3600.D0
        IF(PLTDECSN.EQ.'-')THEN
          DEC0=-DEC0
        ELSEIF(PLTDECSN.EQ.'+')THEN
        ELSE
ccc suponemos que es positiva
ccc          WRITE(*,101)'FATAL ERROR: declination has not sign.'
ccc          STOP
        END IF
        DEC0=DEC0*PI/180.D0
C------------------------------------------------------------------------------
C corregimos de 0.5 pixel de origen (defecto tomado por POSS)
ccc SIN BINNING
ccc        X=(PPO3-XPIXELSZ*DBLE(REAL(CNPIX1-1)+XC+0.5))/1000.D0
ccc        Y=(YPIXELSZ*DBLE(REAL(CNPIX2-1)+YC+0.5)-PPO6)/1000.D0
c tenemos en cuenta el binning
        X=PPO3-XPIXELSZ*DBLE(REAL(CNPIX1-1)+(XC-0.5)*REAL(NBIN)+1.)
        X=X/1000.D0
        Y=YPIXELSZ*DBLE(REAL(CNPIX2-1)+(YC-0.5)*REAL(NBIN)+1.)-PPO6
        Y=Y/1000.D0
C
        XX=AMDX1*X+AMDX2*Y+AMDX3+AMDX4*X*X+AMDX5*X*Y+AMDX6*Y*Y+
     +   AMDX7*(X*X+Y*Y)+AMDX8*X*X*X+AMDX9*X*X*Y+AMDX10*X*Y*Y+
     +   AMDX11*Y*Y*Y+AMDX12*X*(X*X+Y*Y)+AMDX13*X*(X*X+Y*Y)*
     +   (X*X+Y*Y)
        YY=AMDY1*Y+AMDY2*X+AMDY3+AMDY4*Y*Y+AMDY5*X*Y+AMDY6*X*X+
     +   AMDY7*(X*X+Y*Y)+AMDY8*Y*Y*Y+AMDY9*X*Y*Y+AMDY10*X*X*Y+
     +   AMDY11*X*X*X+AMDY12*Y*(X*X+Y*Y)+AMDY13*Y*(X*X+Y*Y)*
     +   (X*X+Y*Y)
        XX=XX/3600.D0                                         !pasamos a grados
        YY=YY/3600.D0                                         !pasamos a grados
        XX=XX*PI/180.D0                                     !pasamos a radianes
        YY=YY*PI/180.D0                                     !pasamos a radianes
        RA1=DATAN((XX/DCOS(DEC0))/(1.D0-YY*DTAN(DEC0)))+RA0
        IF(RA1.LT.0.D0) RA1=RA1+2.D0*PI
        DEC1=DATAN(((YY+DTAN(DEC0))*DCOS(RA1-RA0))
     +   /(1.D0-YY*DTAN(DEC0)))
        RA1=RA1*180.D0/PI                                     !pasamos a grados
        RA1=RA1/15.D0                                          !pasamos a horas
        DEC1=DEC1*180.D0/PI                                   !pasamos a grados
        WRITE(*,100)'(J2000.0) '
        ARHD=DINT(RA1)
        ARMD=(RA1-ARHD)*60.D0
        ARSD=(ARMD-DINT(ARMD))*60.D0
        ARMD=DINT(ARMD)
        IF(DEC1.LT.0.D0)THEN
          DECSIGD='-'
          DEC1=-DEC1
        ELSE
          DECSIGD='+'
        END IF
        DECDD=DINT(DEC1)
        DECMD=(DEC1-DECDD)*60.D0
        DECSD=(DECMD-DINT(DECMD))*60.D0
        DECMD=DINT(DECMD)
        WRITE(*,'(A,I2.2,1X,I2.2,1X,F5.2,A,I2.2,1X,I2.2,1X,F5.2)')
     +   'R.A.: ',INT(ARHD),INT(ARMD),REAL(ARSD),
     +   '  DEC.: '//DECSIGD,INT(DECDD),INT(DECMD),REAL(DECSD)
        IF(DECSIGD.EQ.'-') DEC1=-DEC1
C
        ARH=INT(ARHD)
        ARM=INT(ARMD)
        ARS=REAL(ARSD)
        DECSIG=DECSIGD
        DECD=INT(DECDD)
        DECM=INT(DECMD)
        DECS=REAL(DECSD)
C
        CALL SUBPRECE(2000.D0,RA1,DEC1,1950.D0,RA1,DEC1)
C
        WRITE(*,100)'(B1950.0) '
        ARHD=DINT(RA1)
        ARMD=(RA1-ARHD)*60.D0
        ARSD=(ARMD-DINT(ARMD))*60.D0
        ARMD=DINT(ARMD)
        IF(DEC1.LT.0.D0)THEN
          DECSIGD='-'
          DEC1=-DEC1
        ELSE
          DECSIGD='+'
        END IF
        DECDD=DINT(DEC1)
        DECMD=(DEC1-DECDD)*60.D0
        DECSD=(DECMD-DINT(DECMD))*60.D0
        DECMD=DINT(DECMD)
        WRITE(*,'(A,I2.2,1X,I2.2,1X,F5.2,A,I2.2,1X,I2.2,1X,F5.2)')
     +   'R.A.: ',INT(ARHD),INT(ARMD),REAL(ARSD),
     +   '  DEC.: '//DECSIGD,INT(DECDD),INT(DECMD),
     +   REAL(DECSD)
C------------------------------------------------------------------------------
100     FORMAT(A,$)
        END
C
C******************************************************************************
C******************************************************************************
C
C Introducimos R.A. and DEC. y calcula posicion del cursor
C
        SUBROUTINE COOR_XY(ARH,ARM,ARS,DECSIG,DECD,DECM,DECS,XC,YC)
        IMPLICIT NONE
        INTEGER ARH,ARM
        REAL ARS
        CHARACTER*1 DECSIG
        INTEGER DECD,DECM
        REAL DECS
        REAL XC,YC
C
        DOUBLE PRECISION PI
        PARAMETER(PI=3.14159265358979323846D0)
C
        INTEGER NBIN
        INTEGER ITER
        INTEGER CNPIX1,CNPIX2
        INTEGER PLTRAH,PLTRAM
        INTEGER PLTDECD,PLTDECM
        REAL PRECISION
        DOUBLE PRECISION PPO3,PPO6
        DOUBLE PRECISION XPIXELSZ,YPIXELSZ
        DOUBLE PRECISION AMDX1,AMDX2,AMDX3,AMDX4,AMDX5,AMDX6,AMDX7
        DOUBLE PRECISION AMDX8,AMDX9,AMDX10,AMDX11,AMDX12,AMDX13
        DOUBLE PRECISION AMDY1,AMDY2,AMDY3,AMDY4,AMDY5,AMDY6,AMDY7
        DOUBLE PRECISION AMDY8,AMDY9,AMDY10,AMDY11,AMDY12,AMDY13
        DOUBLE PRECISION XX,YY,ZZ
        DOUBLE PRECISION DETERMINANTE
        DOUBLE PRECISION X0,Y0,XC0,YC0
        DOUBLE PRECISION X1,Y1,XC1,YC1
        DOUBLE PRECISION RA0,DEC0,RA1,DEC1
        DOUBLE PRECISION PLTRAS,PLTDECS
        CHARACTER*1 PLTDECSN
C
        COMMON/BLKNBIN/NBIN
        COMMON/BLKPOSS1/PPO3,PPO6,XPIXELSZ,YPIXELSZ
        COMMON/BLKPOSS2/AMDX1,AMDX2,AMDX3,AMDX4,AMDX5,AMDX6,AMDX7,
     +   AMDX8,AMDX9,AMDX10,AMDX11,AMDX12,AMDX13
        COMMON/BLKPOSS3/AMDY1,AMDY2,AMDY3,AMDY4,AMDY5,AMDY6,AMDY7,
     +   AMDY8,AMDY9,AMDY10,AMDY11,AMDY12,AMDY13
        COMMON/BLKPOSS4/CNPIX1,CNPIX2
        COMMON/BLKPOSS5/PLTRAH,PLTRAM,PLTDECD,PLTDECM
        COMMON/BLKPOSS6/PLTRAS,PLTDECS
        COMMON/BLKPOSS7/PLTDECSN
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Coordenadas del centro de la placa
        RA0=DBLE(PLTRAH)+DBLE(PLTRAM)/60.D0+PLTRAS/3600.D0
        RA0=RA0*15.D0
        RA0=RA0*PI/180.D0
        DEC0=DBLE(PLTDECD)+DBLE(PLTDECM)/60.D0+PLTDECS/3600.D0
        IF(PLTDECSN.EQ.'-')THEN
          DEC0=-DEC0
        ELSEIF(PLTDECSN.EQ.'+')THEN
        ELSE
ccc suponemos que es positiva
ccc          WRITE(*,101)'FATAL ERROR: declination has not sign.'
ccc          STOP
        END IF
        DEC0=DEC0*PI/180.D0
C------------------------------------------------------------------------------
C Coordenadas del objeto
        RA1=DBLE(ARH)+DBLE(ARM)/60.D0+DBLE(ARS)/3600.D0
        RA1=RA1*15.D0*PI/180.D0
        DEC1=DBLE(DECD)+DBLE(DECM)/60.D0+DBLE(DECS)/3600.D0
        IF(DECSIG.EQ.'-')THEN
          DEC1=-DEC1
        ELSEIF(DECSIG.EQ.'+')THEN
        ELSE
ccc suponemos que es positiva
ccc          WRITE(*,101)'FATAL ERROR: declination has not sign.'
ccc          STOP
        END IF
        DEC1=DEC1*PI/180.D0
C------------------------------------------------------------------------------
C Coordenadas plano del cielo
        ZZ=DSIN(DEC0)*DSIN(DEC1)+DCOS(DEC0)*DCOS(DEC1)*DCOS(RA1-RA0)
        XX=DCOS(DEC1)*DSIN(RA1-RA0)
        YY=DCOS(DEC0)*DSIN(DEC1)-DSIN(DEC0)*DCOS(DEC1)*DCOS(RA1-RA0)
        XX=XX/ZZ
        YY=YY/ZZ
        XX=XX*180.D0/PI                                       !pasamos a grados
        YY=YY*180.D0/PI                                       !pasamos a grados
        XX=XX*3600.D0                               !pasamos a segundos de arco
        YY=YY*3600.D0                               !pasamos a segundos de arco
C------------------------------------------------------------------------------
C primera aproximacion (resolvemos el caso lineal)
        DETERMINANTE=AMDX1*AMDY1-AMDX2*AMDY2
        X0=((XX-AMDX3)*AMDY1-AMDX2*(YY-AMDY3))/DETERMINANTE
        Y0=(AMDX1*(YY-AMDY3)-(XX-AMDX3)*AMDY2)/DETERMINANTE
ccc SIN BINNING
ccc        XC0=REAL((PPO3-X0*1000.D0)/XPIXELSZ)-REAL(CNPIX1-1)-0.5
ccc        YC0=REAL((Y0*1000.D0+PPO6)/YPIXELSZ)-REAL(CNPIX2-1)-0.5
ccc con binning
        XC0=REAL((PPO3-X0*1000.D0)/XPIXELSZ)-REAL(CNPIX1-1)-1.0+
     +   0.5*REAL(NBIN)
        XC0=XC0/REAL(NBIN)
        YC0=REAL((Y0*1000.D0+PPO6)/YPIXELSZ)-REAL(CNPIX2-1)-1.0+
     +   0.5*REAL(NBIN)
        YC0=YC0/REAL(NBIN)
ccc        type*,'>>> ',0,xc0,yc0
C------------------------------------------------------------------------------
C Iteramos mediante el proceso de Seidel (ver Demidovich y Maron, p. 174)
        ITER=0
        PRECISION=0.01                                    !en unidades de pixel
10      ITER=ITER+1
        X1=XX-AMDX2*Y0-AMDX3-AMDX4*X0*X0-AMDX5*X0*Y0-AMDX6*Y0*Y0-
     +   AMDX7*(X0*X0+Y0*Y0)-AMDX8*X0*X0*X0-AMDX9*X0*X0*Y0-
     +   AMDX10*X0*Y0*Y0-AMDX11*Y0*Y0*Y0-AMDX12*X0*(X0*X0+Y0*Y0)-
     +   AMDX13*X0*(X0*X0+Y0*Y0)*(X0*X0+Y0*Y0)
        X1=X1/AMDX1
        Y1=YY-AMDY2*X1-AMDY3-AMDY4*Y0*Y0-AMDY5*X1*Y0-AMDY6*X1*X1-
     +   AMDY7*(X1*X1+Y0*Y0)-AMDY8*Y0*Y0*Y0-AMDY9*X1*Y0*Y0-
     +   AMDY10*X1*X1*Y0-AMDY11*X1*X1*X1-AMDY12*Y0*(X1*X1+Y0*Y0)-
     +   AMDY13*Y0*(X1*X1+Y0*Y0)*(X1*X1+Y0*Y0)
        Y1=Y1/AMDY1
C
        XC1=REAL((PPO3-X1*1000.D0)/XPIXELSZ)-REAL(CNPIX1-1)-1.0+
     +   0.5*REAL(NBIN)
        XC1=XC1/REAL(NBIN)
        YC1=REAL((Y1*1000.D0+PPO6)/YPIXELSZ)-REAL(CNPIX2-1)-1.0+
     +   0.5*REAL(NBIN)
        YC1=YC1/REAL(NBIN)
ccc        type*,'>>> ',iter,xc1,yc1
C
        IF((ABS(XC0-XC1).GE.PRECISION).OR.
     +   (ABS(YC0-YC1).GE.PRECISION))THEN
          X0=X1
          Y0=Y1
          XC0=XC1
          YC0=YC1
          GOTO 10
        END IF
        XC=XC1
        YC=YC1
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C Convierte a J2000. La salida se realiza a traves de los mismos parametros
C de entrada
C
        SUBROUTINE PRECESOBJ(EQUINOX,ARH,ARM,ARS,CSIG,DECD,DECM,DECS)
        IMPLICIT NONE
C
        REAL EQUINOX
        INTEGER ARH,ARM
        REAL ARS
        CHARACTER*1 CSIG
        INTEGER DECD,DECM
        REAL DECS
C
        DOUBLE PRECISION PI
        PARAMETER(PI=3.14159265358979323846D0)
C
        DOUBLE PRECISION RA0,DEC0
        DOUBLE PRECISION RA1,DEC1
        DOUBLE PRECISION ARHD,ARMD,ARSD
        DOUBLE PRECISION DECDD,DECMD,DECSD
C------------------------------------------------------------------------------
        RA0=DBLE(ARH)+DBLE(ARM)/60.D0+DBLE(ARS)/3600.D0
        DEC0=DBLE(DECD)+DBLE(DECM)/60.D0+DBLE(DECS)/3600.D0
        IF(CSIG.EQ.'-')THEN
          DEC0=-DEC0
        ELSEIF(CSIG.EQ.'+')THEN
        ELSE
ccc suponemos que es positiva
ccc          WRITE(*,101)'FATAL ERROR: declination has not sign.'
ccc          STOP
        END IF
C
        CALL SUBPRECE(DBLE(EQUINOX),RA0,DEC0,2000.D0,RA1,DEC1)
C
        ARHD=DINT(RA1)
        ARMD=(RA1-ARHD)*60.D0
        ARSD=(ARMD-DINT(ARMD))*60.D0
        ARMD=DINT(ARMD)
        IF(DEC1.LT.0.D0)THEN
          CSIG='-'
          DEC1=-DEC1
        ELSE
          CSIG='+'
        END IF
        DECDD=DINT(DEC1)
        DECMD=(DEC1-DECDD)*60.D0
        DECSD=(DECMD-DINT(DECMD))*60.D0
        DECMD=DINT(DECMD)
        IF(CSIG.EQ.'-') DEC1=-DEC1              !no es necesario pero mejor asi
C
        ARH=INT(ARHD)
        ARM=INT(ARMD)
        ARS=REAL(ARSD)
        DECD=INT(DECDD)
        DECM=INT(DECMD)
        DECS=REAL(DECSD)
C
        END
C
C******************************************************************************
C
        SUBROUTINE PRINTERROR(ISTATUS)
C Print out the FITSIO error messages to the user
        INTEGER ISTATUS
        CHARACTER ERRTEXT*30,ERRMESSAGE*80
C Check if status is OK (no error); if so, simply return
        IF(ISTATUS.LE.0) RETURN
C Get the text string which describes the error
        CALL FTGERR(ISTATUS,ERRTEXT)
        WRITE(*,'(A,$)')'FITSIO Error Status = '
        WRITE(*,*)ISTATUS
        WRITE(*,'(A)')ERRTEXT
C Read and print out all the error messages on the FITSIO stack
        CALL FTGMSG(ERRMESSAGE)
        DO WHILE(ERRMESSAGE.NE.' ')
          WRITE(*,'(A)') ERRMESSAGE
          CALL FTGMSG(ERRMESSAGE)
        END DO
        END
