C------------------------------------------------------------------------------
C Version 30-November-1998                                      file:mcolumns.f
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
C Program: mcolumns
C Classification: miscellany
C Description: Manipulate columns of data in ASCII files.
C
Comment
C
        PROGRAM MCOLUMNS
        IMPLICIT NONE
C
        INCLUDE 'futils.inc'
        INTEGER TRUEBEG,TRUELEN
        INTEGER READI
        INTEGER READILIM
C
        INTEGER NMAX_FILES
        PARAMETER (NMAX_FILES=5)               !numero maximo de ficheros input
        INTEGER NMAX_COLUMNS
        PARAMETER (NMAX_COLUMNS=40)          !numero maximo de columnas/fichero
        INTEGER LMAX
        PARAMETER (LMAX=30) !numero maximo de caracteres mostrados en nombre de
                                                                      !ficheros
C
        INTEGER NFILES,NDELETE,NFILESTATUS
        INTEGER I,J,K,JOUT,ISPACE,L,IS1,IS2
        INTEGER LF(NMAX_FILES)                               !TRUELEN de FILEIN
        INTEGER NLSKIP(NMAX_FILES)  !numero de lineas a ignorar en cada fichero
        INTEGER NCOL(NMAX_FILES)            !numero de columnas en cada fichero
        INTEGER NROW(NMAX_FILES)               !numero de filas en cada fichero
        INTEGER LFNEW,NLSKIPNEW,NCOLNEW,NROWNEW,NROWMIN,NROWOUT
        INTEGER NOFCOL       !funcion para calcular numero de columnas en LINEA
        INTEGER NCOLOUT,NBLANK
        INTEGER NNCOL(NMAX_COLUMNS)
        INTEGER NNFILE(NMAX_COLUMNS)
        INTEGER NNSEP(NMAX_COLUMNS)
        INTEGER WMIN(NMAX_COLUMNS,NMAX_FILES)         !anchura minima del campo
        INTEGER WMAX(NMAX_COLUMNS,NMAX_FILES)         !anchura maxima del campo
        INTEGER WLMIN(NMAX_FILES)                   !anchura minima de la linea
        INTEGER WLMAX(NMAX_FILES)                   !anchura maxima de la linea
        INTEGER WL1,WL2,WL0
        INTEGER WPREFIX(NMAX_COLUMNS)      !anchura maxima del prefijo (if any)
        CHARACTER*1 COPC,CCONT,CBLANK
        CHARACTER*1 CJUST(NMAX_COLUMNS)         !justificacion del texto: l,c,r
        CHARACTER*1 CJUST_TEMPORAL
        CHARACTER*75 FILEIN(NMAX_FILES)      !nombre de los ficheros de entrada
        CHARACTER*75 NEWFILEIN,OUTFILE
        CHARACTER*10 PREFIX(NMAX_COLUMNS)
        CHARACTER*255 BLANKLINE,BLANKFIELD,CDUMMY
        CHARACTER*255 LINEA                       !linea a leer en cada fichero
        CHARACTER*255 LINEAIN(NMAX_FILES)         !linea a leer en cada fichero
        CHARACTER*255 DATO(NMAX_COLUMNS,NMAX_FILES)     !dato (columna,fichero)
        LOGICAL LLOOP,LOGFILE,LROWOK,LANYROW
        LOGICAL IFCOLOUT(NMAX_COLUMNS)   !.TRUE. si la columna ha sido definida
        CHARACTER*255 DATOLINEA(NMAX_COLUMNS)
C
        COMMON/BLKDATOLINEA/DATOLINEA
        COMMON/BLKNFILES/NFILES
C------------------------------------------------------------------------------
        NBLANK=0 !evitamos un warning de compilacion
C Seleccionamos los ficheros a utilizar
        NFILES=0
        LLOOP=.TRUE.
        CBLANK='%'
        DO I=1,NMAX_COLUMNS
          CJUST(I)='l'
        END DO
        DO I=1,LEN(BLANKLINE)
          BLANKLINE(I:I)=' '
        END DO
C
        DO I=1,NMAX_COLUMNS
          IFCOLOUT(I)=.FALSE.
        END DO
C------------------------------------------------------------------------------
C listamos ficheros input
10      WRITE(*,100)'[1;1f[J'
        WRITE(*,100)'[32m'
        WRITE(*,'(A,I1)') '* List of input ASCII files: ',NFILES
        WRITE(*,100)'[33m'
        IF(NFILES.GT.0)THEN
          DO I=1,NFILES
            IF(LF(I).LT.LMAX)THEN
              WRITE(*,'(A,I1,A,A,A,A,$)') '#',I,'> "',
     +         FILEIN(I)(1:LF(I)),'"',BLANKLINE(1:LMAX-LF(I))
            ELSE
              WRITE(*,'(A,I1,A,$)') '#',I,'> "'
              WRITE(*,100)'[31m'
              WRITE(*,100) '...'
              WRITE(*,100)'[33m'
              WRITE(*,'(A,A,$)') FILEIN(I)(LF(I)-LMAX+4:LF(I)),'"'
            END IF
            WRITE(*,100)' '
            WRITE(*,100)'[7m'
            WRITE(*,'(A,I4.4,$)')'skip: ',NLSKIP(I)
            WRITE(*,100)'[0m'
            WRITE(*,100)'[33m'
            WRITE(*,100)' '
            WRITE(*,100)'[7m'
            WRITE(*,'(A,I4.4,$)')'ncol: ',NCOL(I)
            WRITE(*,100)'[0m'
            WRITE(*,100)'[33m'
            WRITE(*,100)' '
            WRITE(*,100)'[7m'
            WRITE(*,'(A,I4.4,$)')'nrow: ',NROW(I)
            WRITE(*,100)'[0m'
            WRITE(*,100)'[33m'
            WRITE(*,100)' '
            WRITE(*,100)'[7m'
            WRITE(*,'(A,I3.3,A,I3.3)')'w:',WLMIN(I),'/',WLMAX(I)
            WRITE(*,100)'[0m'
            WRITE(*,100)'[33m'
          END DO
        END IF
C
        IF(NFILES.LT.NMAX_FILES)THEN
          DO I=NFILES+1,NMAX_FILES
            WRITE(*,100)'[33m'
            WRITE(*,'(A,I2.2,A)') '#',I,'> ***file not defined***'
          END DO
        END IF
C------------------------------------------------------------------------------
C listamos status de columnas en fichero output
        WRITE(*,*)
        WRITE(*,100)'[32m'
        WRITE(*,101)'* Status of columns in output ASCII file:'
        WRITE(*,100)'[35m'
        DO I=1,NMAX_COLUMNS
          CALL LOCATE(I,1)
          IF(IFCOLOUT(I))THEN
            WRITE(*,100)'[32m'
            WRITE(*,100)'[7m'
            IF(NNFILE(I).EQ.0)THEN
              IF(TRUELEN(PREFIX(I)).GT.0)THEN
                WRITE(*,'(A,I2.2,A,I1,A,I1,A,I2.2,A)')'#',I,
     +           '>F',NNFILE(I),'+p',WPREFIX(I),' S',NNSEP(I),CJUST(I)
              ELSE
                WRITE(*,'(A,I2.2,A,I1,A,I1,A,I2.2,A)')'#',I,
     +           '>F',NNFILE(I),'-p',WPREFIX(I),' S',NNSEP(I),CJUST(I)
              END IF
            ELSE
              WRITE(*,'(A,I2.2,A,I1,A,I2.2,A,I2.2,A)')'#',I,
     +         '>F',NNFILE(I),' C',NNCOL(I),' S',NNSEP(I),CJUST(I)
            END IF
            WRITE(*,100)'[0m'
          ELSE
            WRITE(*,100)'[35m'
            WRITE(*,'(A,I2.2)')'#',I
          END IF
        END DO
C------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100)'[36m'
        WRITE(*,100)'[i]'
        WRITE(*,100)'[0m'
        WRITE(*,100)'nclude new file, '
        WRITE(*,100)'[36m'
        WRITE(*,100)'[d]'
        WRITE(*,100)'[0m'
        WRITE(*,100)'elete file, '
        WRITE(*,100)'[36m'
        WRITE(*,100)'[c]'
        WRITE(*,100)'[0m'
        WRITE(*,100)'reate output file, '
        WRITE(*,100)'[36m'
        WRITE(*,100)'[s]'
        WRITE(*,100)'[0m'
        WRITE(*,101)'tatus input file, '
        WRITE(*,100)'[36m'
        WRITE(*,100)'[+]'
        WRITE(*,100)'[0m'
        WRITE(*,100)'add output column, '
        WRITE(*,100)'[36m'
        WRITE(*,100)'[-]'
        WRITE(*,100)'[0m'
        WRITE(*,100)'delete output column, '
        WRITE(*,100)'[36m'
        WRITE(*,100)'[q]'
        WRITE(*,100)'[0m'
        WRITE(*,101)'uit'
        WRITE(*,100)'Option (i/d/c/s/+/-/q) '
        WRITE(*,100)'[0m'
        COPC(1:1)=READC('@','iIdDcCsS+-qQ')
        CALL CHLOWER(COPC)
        WRITE(*,100)'[0m'
C------------------------------------------------------------------------------
C quit
        IF(COPC.EQ.'q')THEN
          STOP
        ELSEIF((COPC.EQ.'d').OR.(COPC.EQ.'s'))THEN
          IF(NFILES.EQ.0)THEN
            WRITE(*,100)'[31m'
            WRITE(*,100)'ERROR: No. of files = 0. '
            WRITE(*,101)'You must define at lest one input file!'
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
            WRITE(*,100)'[0m'
            GOTO 10
          END IF
        END IF
C------------------------------------------------------------------------------
C include file
        IF(COPC.EQ.'i')THEN
          IF(NFILES.EQ.NMAX_FILES)THEN
            WRITE(*,100)'[31m'
            WRITE(*,100)'* WARNING: No. of files = maximum number. '
            WRITE(*,101)'You must delete one file before'
            WRITE(*,101)'inserting a new one.'
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
            WRITE(*,100)'[0m'
          ELSE
            WRITE(*,100)'New input ASCII file name? '
            READ(*,101)NEWFILEIN
            LFNEW=TRUELEN(NEWFILEIN)
            IF(LFNEW.EQ.0) GOTO 10
            INQUIRE(FILE=NEWFILEIN,EXIST=LOGFILE)
            IF(LOGFILE)THEN
              IF(NFILES.GT.0)THEN          !chequeamos que no repetimos fichero
                DO I=1,NFILES
                  IF(NEWFILEIN(1:LFNEW).EQ.FILEIN(I)(1:LF(I)))THEN
                    WRITE(*,100)'[31m'
                    WRITE(*,101)'ERROR: this file is in the list.'
                    WRITE(*,100)'Press <CR> to continue...'
                    READ(*,*)
                    WRITE(*,100)'[0m'
                    GOTO 10
                  END IF
                END DO
              END IF
              WRITE(*,100)'No. of lines to be skipped '
              NLSKIPNEW=READILIM('0',0,9999)
C comprobamos que el fichero tiene lineas suficientes
              OPEN(10,FILE=NEWFILEIN,STATUS='OLD',FORM='FORMATTED')
              IF(NLSKIPNEW.GT.0)THEN
                DO I=1,NLSKIPNEW
                  READ(10,*,END=20)
                END DO
              END IF
              GOTO 21
20            CLOSE(10)
              WRITE(*,100)'[31m'
              WRITE(*,100)'ERROR: this file does not contain enough '
              WRITE(*,101)'lines to be skipped!'
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
              WRITE(*,100)'[0m'
              GOTO 10
21            READ(10,101,END=22) LINEA
              GOTO 23
22            CLOSE(10)
              WRITE(*,100)'[31m'
              WRITE(*,100)'ERROR: this file does not contain data '
              WRITE(*,101)'lines!'
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
              WRITE(*,100)'[0m'
              GOTO 10
C averiguamos el numero de columnas
23            WL1=TRUEBEG(LINEA)
              WL2=TRUELEN(LINEA)
              WLMIN(NFILES+1)=WL2-WL1+1
              WLMAX(NFILES+1)=WL2-WL1+1
              NCOLNEW=NOFCOL(LINEA)
              IF(NCOLNEW.EQ.0)THEN
                CLOSE(10)
                WRITE(*,100)'[31m'
                WRITE(*,101)'ERROR: empty data line!'
                WRITE(*,100)'Press <CR> to continue...'
                READ(*,*)
                WRITE(*,100)'[0m'
                GOTO 10
              ELSEIF(NCOLNEW.GT.NMAX_COLUMNS)THEN
                CLOSE(10)
                WRITE(*,100)'[31m'
                WRITE(*,100)'ERROR: no. of data columns is too large '
                WRITE(*,'(A,I2.2,A)')'(NCOL > ',NMAX_COLUMNS,')'
                WRITE(*,100)'Press <CR> to continue...'
                READ(*,*)
                WRITE(*,100)'[0m'
                GOTO 10
              END IF
C averiguamos el numero de lineas con datos (ya hemos leido una!)
              NROWNEW=1
              DO J=1,NCOLNEW
                L=TRUELEN(DATOLINEA(J))
                WMIN(J,NFILES+1)=L !determinaremos la longitud minima del campo
                WMAX(J,NFILES+1)=L !determinaremos la longitud maxima del campo
              END DO
              WL1=TRUEBEG(LINEA)
              WL2=TRUELEN(LINEA)
              WL0=WL2-WL1+1
              IF(WL0.LT.WLMIN(NFILES+1)) WLMIN(NFILES+1)=WL0
              IF(WL0.GT.WLMAX(NFILES+1)) WLMAX(NFILES+1)=WL0
24            READ(10,101,END=25) LINEA
              NROWNEW=NROWNEW+1
              IF(NOFCOL(LINEA).NE.NCOLNEW)THEN
                CLOSE(10)
                WRITE(*,100)'[31m'
                WRITE(*,100)'ERROR: no. of data columns changes in '
                WRITE(*,'(A,I2.2)')'data line number: ',NROWNEW
                WRITE(*,100)'Expected number of columns: '
                WRITE(*,*) NCOLNEW
                WRITE(*,100)'Actual number of columns..: '
                WRITE(*,*) NOFCOL(LINEA)
                WRITE(*,100)'Press <CR> to continue...'
                READ(*,*)
                WRITE(*,100)'[0m'
                GOTO 10
              END IF
              DO J=1,NCOLNEW
                L=TRUELEN(DATOLINEA(J))
                IF(L.LT.WMIN(J,NFILES+1)) WMIN(J,NFILES+1)=L
                IF(L.GT.WMAX(J,NFILES+1)) WMAX(J,NFILES+1)=L
              END DO
              GOTO 24
25            CLOSE(10)
C el fichero parece ser correcto: lo almacenamos en la lista de ficheros
              NFILES=NFILES+1
              FILEIN(NFILES)=NEWFILEIN
              LF(NFILES)=TRUELEN(NEWFILEIN)
              NLSKIP(NFILES)=NLSKIPNEW
              NCOL(NFILES)=NCOLNEW
              NROW(NFILES)=NROWNEW
            ELSE
              WRITE(*,100)'[31m'
              WRITE(*,101)'ERROR: this file does not exist.'
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
              WRITE(*,100)'[0m'
            END IF
          END IF
C------------------------------------------------------------------------------
C delete file
        ELSEIF(COPC.EQ.'d')THEN
          WRITE(*,100)'No. of file to be removed (0=none) '
          NDELETE=READILIM('@',0,NFILES)
          IF(NDELETE.EQ.0)THEN
            GOTO 10
          ELSEIF(NDELETE.EQ.NFILES)THEN
            NFILES=NFILES-1
          ELSE
            DO I=NDELETE,NFILES-1
              FILEIN(I)=FILEIN(I+1)
              LF(I)=LF(I+1)
              NLSKIP(I)=NLSKIP(I+1)
              NCOL(I)=NCOL(I+1)
              NROW(I)=NROW(I+1)
              DO J=1,NCOL(I+1)
                WMIN(J,I)=WMIN(J,I+1)
                WMAX(J,I)=WMAX(J,I+1)
              END DO
            END DO
            NFILES=NFILES-1
          END IF
C------------------------------------------------------------------------------
C create output file
        ELSEIF(COPC.EQ.'c')THEN
C comprobamos si hay definida alguna columna de salida
          LANYROW=.FALSE.
          DO I=1,NMAX_COLUMNS
            IF(IFCOLOUT(I)) LANYROW=.TRUE.
          END DO
          IF(LANYROW)THEN
          ELSE
            WRITE(*,100)'[31m'
            WRITE(*,101)'ERROR: no. of defined output columns=0!'
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
            WRITE(*,100)'[0m'
            GOTO 10
          END IF
C comprobamos longitud de datos en ficheros
          LROWOK=.TRUE.
          NROWMIN=NROW(1)
          IF(NFILES.GT.1)THEN
            DO I=2,NFILES
              IF(NROW(I).NE.NROWMIN)THEN
                LROWOK=.FALSE.
                IF(NROW(I).LT.NROWMIN) NROWMIN=NROW(I)
              END IF
            END DO
          END IF
          IF(LROWOK)THEN
            NROWOUT=NROWMIN
          ELSE
            WRITE(*,100)'[31m'
            WRITE(*,101)'No. of rows in input files are different.'
            WRITE(*,100)'Continue anyway (y/n) '
            CCONT(1:1)=READC('y','yn')
            WRITE(*,100)'[0m'
            IF(CCONT.EQ.'n') GOTO 10
            WRITE(*,100)'[22;1f'
            WRITE(*,100)'[K'
            WRITE(*,100)'[21;1f'
            WRITE(*,100)'[K'
            WRITE(*,100)'No. of rows in output file '
            WRITE(CDUMMY,*) NROWMIN
            NROWOUT=READI(CDUMMY)
            IF(NROWOUT.LT.1)THEN
              WRITE(*,100)'[31m'
              WRITE(*,101)'ERROR: no. of rows < 1!'
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
              WRITE(*,100)'[0m'
              GOTO 10
            END IF
            IF(NROWOUT.GT.NROWMIN)THEN
              WRITE(*,100)'Blank field character to extend columns '
              CBLANK(1:1)=READC(CBLANK,'@')
              WRITE(*,100)'[1A'
              WRITE(*,100)'No. of blank field characters in '//
     +         'each column '
              NBLANK=READILIM('255',1,255)
            END IF
          END IF
C pedimos nombre de fichero de salida
          WRITE(*,100)'Output file name'
          OUTFILE(1:75)=READC('@','@')
          INQUIRE(FILE=OUTFILE,EXIST=LOGFILE)
          IF(LOGFILE)THEN
            WRITE(*,100)'[31m'
            WRITE(*,100)'ERROR: this file already exist. You can '
            WRITE(*,101)'not overwrite it.'
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
            WRITE(*,100)'[0m'
            GOTO 10
          END IF
C generamos fichero de salida
          OPEN(20,FILE=OUTFILE,STATUS='NEW',FORM='FORMATTED')
          DO I=1,NFILES
            OPEN(20+I,FILE=FILEIN(I),STATUS='OLD',FORM='FORMATTED')
            IF(NLSKIP(I).GT.0)THEN         !en caso necesario, ignoramos lineas
              DO K=1,NLSKIP(I)
                READ(20+I,*)
              END DO
            END IF
          END DO
C
          DO I=1,LEN(BLANKFIELD)    !definimos la linea BLANK por si hace falta
            BLANKFIELD(I:I)=CBLANK
          END DO
C
          DO K=1,NROWOUT                !para cada fila en el fichero de salida
            DO I=1,NFILES       !leemos informacion en todos los ficheros input
              WL1=1
              IF(K.GT.NROW(I))THEN              !si tenemos que alargar columna
                DO J=1,NCOL(I)
                  IF(NBLANK.LE.WMAX(J,I))THEN
                    DATO(J,I)=BLANKFIELD(1:NBLANK)
                  ELSE
                    DATO(J,I)=BLANKFIELD(1:WMAX(J,I))
                  END IF
                  WL0=TRUELEN(DATO(J,I))
                  WL2=WL1+WL0-1
                  LINEAIN(I)(WL1:WL2)=DATO(J,I)(1:WL0)
                  WL1=WL2+2                       !dejamos 2 espacios en blanco
                END DO
              ELSE
                READ(20+I,101)LINEA       !leemos una linea del fichero input I
                NCOLNEW=NOFCOL(LINEA)    !numero de columnas en fichero input I
                DO J=1,NCOLNEW  !guardamos la informacion de todas las columnas
                  DATO(J,I)=DATOLINEA(J) !(esta informacion viene de un COMMON)
                END DO
                LINEAIN(I)=LINEA
              END IF
            END DO
            DO JOUT=1,NMAX_COLUMNS   !recorremos todas las columnas del OUTFILE
              IF(IFCOLOUT(JOUT))THEN           !si ha sido definida, escribimos
                IF(NNFILE(JOUT).EQ.0)THEN        !insertamos numero correlativo
                  WL1=TRUEBEG(PREFIX(JOUT))
                  WL2=TRUELEN(PREFIX(JOUT))
                  WRITE(CDUMMY,'(I9.9)')K
                  WRITE(20,100)PREFIX(JOUT)(WL1:WL2)
                  IF(WPREFIX(JOUT).GT.0)THEN
                    WRITE(20,100)CDUMMY(9-WPREFIX(JOUT)+1:9)
                  END IF
                  IF(NNSEP(JOUT).GT.0)THEN
                    DO ISPACE=1,NNSEP(JOUT)
                      WRITE(20,100)' '
                    END DO
                  END IF
                ELSE
                  IF(NNCOL(JOUT).EQ.0)THEN          !insertamos todo el fichero
                    WL1=TRUEBEG(LINEAIN(NNFILE(JOUT)))
                    WL2=TRUELEN(LINEAIN(NNFILE(JOUT)))
                    WL0=WL2-WL1+1
                    IF(WL0.LT.WLMAX(NNFILE(JOUT)))THEN
                      IF(CJUST(JOUT).EQ.'l')THEN
                        IS1=0
                        IS2=WLMAX(NNFILE(JOUT))-WL0
                      ELSEIF(CJUST(JOUT).EQ.'c')THEN
                        IS1=(WLMAX(NNFILE(JOUT))-WL0)/2
                        IS2=(WLMAX(NNFILE(JOUT))-WL0)-IS1
                      ELSE
                        IS1=WLMAX(NNFILE(JOUT))-WL0
                        IS2=0
                      END IF
                    ELSE
                      IS1=0
                      IS2=0
                    END IF
                    IF(IS1.GT.0)THEN
                      DO ISPACE=1,IS1
                        WRITE(20,100)' '
                      END DO
                    END IF
                    WRITE(20,100)LINEAIN(NNFILE(JOUT))(WL1:WL2)
                    IF(IS2.GT.0)THEN
                      DO ISPACE=1,IS2
                        WRITE(20,100)' '
                      END DO
                    END IF
                    IF(NNSEP(JOUT).GT.0)THEN
                      DO ISPACE=1,NNSEP(JOUT)
                        WRITE(20,100)' '
                      END DO
                    END IF
                  ELSE                     !insertamos solo la columna indicada
                    L=TRUELEN(DATO(NNCOL(JOUT),NNFILE(JOUT)))
                    IF(L.LT.WMAX(NNCOL(JOUT),NNFILE(JOUT)))THEN     !rellenamos
                      IF(CJUST(JOUT).EQ.'l')THEN
                        IS1=0
                        IS2=WMAX(NNCOL(JOUT),NNFILE(JOUT))-L
                      ELSEIF(CJUST(JOUT).EQ.'c')THEN
                        IS1=(WMAX(NNCOL(JOUT),NNFILE(JOUT))-L)/2
                        IS2=(WMAX(NNCOL(JOUT),NNFILE(JOUT))-L)-IS1
                      ELSE
                        IS1=WMAX(NNCOL(JOUT),NNFILE(JOUT))-L
                        IS2=0
                      END IF
                    ELSE
                      IS1=0
                      IS2=0
                    END IF
                    IF(IS1.GT.0)THEN
                      DO ISPACE=1,IS1
                        WRITE(20,100)' '
                      END DO
                    END IF
                    WRITE(20,100)DATO(NNCOL(JOUT),NNFILE(JOUT))(1:L)
                    IF(IS2.GT.0)THEN
                      DO ISPACE=1,IS2
                        WRITE(20,100)' '
                      END DO
                    END IF
                    IF(NNSEP(JOUT).GT.0)THEN
                      DO ISPACE=1,NNSEP(JOUT)
                        WRITE(20,100)' '
                      END DO
                    END IF
                  END IF
                END IF
              END IF
            END DO
            WRITE(20,*)                                        !cambio de linea
          END DO
C
          DO I=0,NFILES
            CLOSE(20+I)
          END DO
C------------------------------------------------------------------------------
C show status of input file
        ELSEIF(COPC.EQ.'s')THEN
          WRITE(*,100)'No. of input file to show status '
          NFILESTATUS=READILIM('@',1,NFILES)
          CALL LOCATE(1,1)
          WRITE(*,100)'[J'
          WRITE(*,100)'[1A'
          WRITE(*,100)'[32m'
          WRITE(*,'(A,I1)')'* Status of columns in input ASCII file #',
     +     NFILESTATUS
          WRITE(*,100)'[35m'
          DO I=1,NMAX_COLUMNS
            CALL LOCATE(I,1)
            IF(I.LE.NCOL(NFILESTATUS))THEN
              WRITE(*,100)'[32m'
              WRITE(*,100)'[7m'
              WRITE(*,'(A,I2.2,A,I3.3,A,I3.3)')'#',I,
     +         '> W:',WMIN(I,NFILESTATUS),'/',WMAX(I,NFILESTATUS)
              WRITE(*,100)'[0m'
            ELSE
              WRITE(*,100)'[35m'
              WRITE(*,'(A,I2.2)')'#',I
              WRITE(*,100)'[0m'
            END IF
          END DO
          WRITE(*,100)'[36m'
          WRITE(*,100)'Press <CR> to continue...'
          READ(*,*)
          WRITE(*,100)'[0m'
C------------------------------------------------------------------------------
C define column in output file
        ELSEIF(COPC.EQ.'+')THEN
          WRITE(*,100)'[35m'   !color
C
          WRITE(*,100)'[21;1f' !mueve cursor a inicio de linea de abajo
          WRITE(*,100)'[K'     !borra hasta fin de linea
          WRITE(*,100)'Number of column in output file '
          NCOLOUT=READILIM('@',1,NMAX_COLUMNS)
          IFCOLOUT(NCOLOUT)=.TRUE.
          CALL LOCATE(NCOLOUT,1)
          WRITE(*,100)'[7m'   !color
          WRITE(*,'(A,I2.2,A,$)')'#',NCOLOUT,'>'
          WRITE(*,100)'[0m'   !color
          WRITE(*,100)'[35m'   !color
C
          WRITE(*,100)'[21;1f' !mueve cursor a inicio de linea de abajo
          WRITE(*,100)'[K'     !borra hasta fin de linea
          WRITE(*,100)'File number (0=insert line number) '
          NNFILE(NCOLOUT)=READILIM('@',0,NFILES)
          CALL LOCATE(NCOLOUT,5)
          WRITE(*,100)'[7m'   !color
          WRITE(*,'(A,I1)')'F',NNFILE(NCOLOUT)
          WRITE(*,100)'[0m'   !color
          WRITE(*,100)'[35m'   !color
C
          IF(NNFILE(NCOLOUT).EQ.0)THEN                     !pedimos nombre raiz
            WRITE(*,100)'[21;1f' !linea de abajo
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'Prefix name (maximum 10 characters)? '
            READ(*,101) PREFIX(NCOLOUT)
            WRITE(*,100)'[21;1f' !linea de abajo
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'No. of digits '
            WPREFIX(NCOLOUT)=READILIM('@',0,6)
            CALL LOCATE(NCOLOUT,8)
            WRITE(*,100)'[7m'   !color
            IF(TRUELEN(PREFIX(NCOLOUT)).GT.0)THEN
              WRITE(*,'(A,I1)')'+p',WPREFIX(NCOLOUT)
            ELSE
              WRITE(*,'(A,I1)')'-p',WPREFIX(NCOLOUT)
            END IF
            WRITE(*,100)'[0m'   !color
            WRITE(*,100)'[35m'   !color
            NNCOL(NCOLOUT)=0
          ELSE                        !insertamos columna de fichero solicitado
            WRITE(*,100)'[21;1f' !linea de abajo
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'Column number in input file (0=all) '
            NNCOL(NCOLOUT)=READILIM('@',0,NCOL(NNFILE(NCOLOUT)))
            CALL LOCATE(NCOLOUT,8)
            WRITE(*,100)'[7m'   !color
            WRITE(*,'(A,I2.2)')'C',NNCOL(NCOLOUT)
            WRITE(*,100)'[0m'   !color
            WRITE(*,100)'[35m'   !color
            PREFIX(NCOLOUT)='UNDEFINED'
            WPREFIX(NCOLOUT)=0
          END IF
C
          WRITE(*,100)'[21;1f' !linea de abajo
          WRITE(*,100)'[K'     !borra hasta fin de linea
          WRITE(*,100)'Column separator (number of blank spaces) '
          NNSEP(NCOLOUT)=READILIM('2',1,20)
          CALL LOCATE(NCOLOUT,12)
          WRITE(*,100)'[7m'   !color
          WRITE(*,'(A,I2.2)')'S',NNSEP(NCOLOUT)
          WRITE(*,100)'[0m'   !color
          WRITE(*,100)'[35m'   !color
C
          IF(NNFILE(NCOLOUT).EQ.0)THEN                     !pedimos nombre raiz
            CJUST(NCOLOUT)='n'
          ELSE
            WRITE(*,100)'[21;1f' !linea de abajo
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'Field justification [l]eft, [c]enter'//
     +       ' or [r]ight '
            CJUST_TEMPORAL(1:1)=READC(CJUST(NCOLOUT),'lcr')
            CJUST(NCOLOUT)=CJUST_TEMPORAL
            CALL LOCATE(NCOLOUT,8)
            WRITE(*,100)'[7m'   !color
            WRITE(*,'(A,I2.2)')'C',NNCOL(NCOLOUT)
            WRITE(*,100)'[0m'   !color
            WRITE(*,100)'[35m'   !color
          END IF
C
          WRITE(*,100)'[21;1f' !linea de abajo
          WRITE(*,100)'[K'     !borra hasta fin de linea
          CALL LOCATE(1,1)
C------------------------------------------------------------------------------
C delete column in output file
        ELSEIF(COPC.EQ.'-')THEN
          WRITE(*,100)'[35m'   !color
C
          WRITE(*,100)'[21;1f' !mueve cursor a inicio de linea de abajo
          WRITE(*,100)'[K'     !borra hasta fin de linea
          WRITE(*,100)'Number of column to be deleted '
          NCOLOUT=READILIM('@',1,NMAX_COLUMNS)
          IF(IFCOLOUT(NCOLOUT))THEN
            IFCOLOUT(NCOLOUT)=.FALSE.
          ELSE
            WRITE(*,100)'[31m'
            WRITE(*,100)'[21;1f'
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'ERROR: this column has not been defined!'
            WRITE(*,100)'[22;1f'
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'Press <CR> to continue...'
            READ(*,*)
            WRITE(*,100)'[21;1f'
            WRITE(*,100)'[K'     !borra hasta fin de linea
            WRITE(*,100)'[22;1f'
            WRITE(*,100)'[K'     !borra hasta fin de linea
          END IF
C
          WRITE(*,100)'[21;1f' !linea de abajo
          WRITE(*,100)'[K'     !borra hasta fin de linea
          CALL LOCATE(1,1)
C------------------------------------------------------------------------------
        END IF
        GOTO 10
C------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        INTEGER FUNCTION NOFCOL(LINEA)
        IMPLICIT NONE
        CHARACTER*255 LINEA
C
        INTEGER TRUEBEG,TRUELEN
C
        INTEGER NMAX_COLUMNS
        PARAMETER (NMAX_COLUMNS=40)          !numero maximo de columnas/fichero
C
        INTEGER L1,L2
        INTEGER NCOL,NEXT
        CHARACTER*255 RESTO
        CHARACTER*255 DATOLINEA(NMAX_COLUMNS)
C
        COMMON/BLKDATOLINEA/DATOLINEA
C------------------------------------------------------------------------------
C caso trivial
        IF(TRUELEN(LINEA).EQ.0)THEN
          NOFCOL=0
          RETURN
        END IF
C
        RESTO=LINEA     !trabajamos con la cadena RESTO para no modificar LINEA
C------------------------------------------------------------------------------
        NCOL=1            !almacenaremos en esta variable el numero de columnas
        L1=TRUEBEG(RESTO)                  !primer elemento valido de la cadena
        L2=TRUELEN(RESTO)                  !ultimo caracter valido de la cadena
C
10      NEXT=INDEX(RESTO(L1:L2),' ')               !siguiente espacio en blanco
        IF(NEXT.EQ.0)THEN        !ya no hay mas espacios en blanco y regresamos
          DATOLINEA(NCOL)=RESTO(L1:L2)
          NOFCOL=NCOL
          RETURN
        END IF
        DATOLINEA(NCOL)=RESTO(L1:L1+NEXT-2)
        NCOL=NCOL+1                 !numero de la siguiente columna a encontrar
        L1=L1+NEXT
        L1=L1+TRUEBEG(RESTO(L1:L2))-1
        GOTO 10
C------------------------------------------------------------------------------
        END
C
C******************************************************************************
C Coloca el cursor en la posicion adecuada para la columna NCOL, con un offset
C hacia la derecha de OFFSET caracteres
        SUBROUTINE LOCATE(NCOL,OFFSET)
        IMPLICIT NONE
        INTEGER NCOL,OFFSET
C
        INTEGER NROWTEXT,NCOLTEXT
        INTEGER NFILES
C
        COMMON/BLKNFILES/NFILES
C------------------------------------------------------------------------------
        NROWTEXT=1+(NCOL-1)/5
        NCOLTEXT=NCOL-(NROWTEXT-1)*5
        WRITE(*,'(A,I2.2,A,I2.2,A,$)')
     +   '[',8+NROWTEXT,";",(NCOLTEXT-1)*16+OFFSET,'f'
C
        END
