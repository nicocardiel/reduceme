C------------------------------------------------------------------------------
C Version 16-April-1998                                          file:adnhand.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: adnhand
C Classification: arithmetic & manipulations
C Description: Create a new spectrum by adding spectra from different images 
C (with interactive monitoring of the S/N ratio).
C
Comment
C
        PROGRAM ADNHAND
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        REAL READF
C
        INTEGER I,J,L
        INTEGER N,N1,N2,N0
        INTEGER NCHAN0,NC1,NC2
        INTEGER NTERMS,NFIT
        INTEGER ITERM,NTERM,IDN(MAX_ID_RED)
        REAL A(NCMAX,NSMAX),S(NCMAX),FS(NCMAX),FS0(NCMAX)
        REAL ERR(NCMAX,NSMAX),ES(NCMAX),FES(NCMAX),FES0(NCMAX)
        REAL FACTOR,XC,YC
        REAL Y(NSMAX),SY(NSMAX),ESY(NSMAX)
        REAL XMIN,XMAX,YMIN,YMAX,YMIN0,YMAX0,DX,DY,X(NCMAX)
        REAL SPD(NCMAX),XFIT(NCMAX),YFIT(NCMAX),COEFF(20),CHISQR
        DOUBLE PRECISION SIGMA,MEAN
        CHARACTER*1 CERR,CNEW,CADD,CEXPAND,CCONF,CSAVE,CMOUSE,CH
        CHARACTER*15 CMEAN1,CMEAN2,CMEAN3
        CHARACTER*15 CSNRA1,CSNRA2,CSNRA3
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
        LOGICAL IFSCAN(NSMAX),LFIRSTIMAGE
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKNC/NC1,NC2
C------------------------------------------------------------------------------
        THISPROGRAM='adnhand'
        CALL WELCOME('16-April-1998')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('y','yn')
C
        DO J=1,NCMAX
          X(J)=REAL(J)
        END DO
        DO J=1,NCMAX
          FS(J)=0.
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCMAX
            FES(J)=0.
          END DO
        END IF
C
        NCHAN0=0
        LFIRSTIMAGE=.TRUE.
C------------------------------------------------------------------------------
C Pedimos salida grafica
        CALL PIDEGTER(NTERM,IDN,LCOLOR)
C------------------------------------------------------------------------------
C Introducimos fichero a utilizar
10      WRITE(*,*)
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        IF(NCHAN0.EQ.0) NCHAN0=NCHAN
        IF(NCHAN.NE.NCHAN0)THEN
          WRITE(*,101)'ERROR: invalid NCHAN value in last image.'
          CLOSE(20)
          GOTO 80
        END IF
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.)!....match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C------------------------------------------------------------------------------
C Examinamos corte espacial
        DO I=1,NSCAN
          Y(I)=REAL(I)
        END DO
        DO I=1,NSCAN
          SY(I)=0.
        END DO
        DO I=1,NSCAN
          ESY(I)=0.
        END DO
C
        DO I=1,NSCAN
          DO J=1,NCHAN
            SY(I)=SY(I)+A(J,I)
          END DO
          SY(I)=SY(I)/REAL(NCHAN)
        END DO
C
        IF(CERR.EQ.'y')THEN
          DO I=1,NSCAN
            DO J=1,NCHAN
              ESY(I)=ESY(I)+ERR(J,I)
            END DO
            ESY(I)=ESY(I)/REAL(NCHAN)
          END DO
        END IF
C
        XMIN=1.
        XMAX=REAL(NSCAN)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
C
        CALL FINDMM(NSCAN,SY,YMIN,YMAX)
        CALL FINDMM(NSCAN,ESY,YMIN0,YMAX0)
        IF(YMIN0.LT.YMIN) YMIN=YMIN0
        IF(YMAX0.GT.YMAX) YMAX=YMAX0
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NSCAN,Y,SY,.TRUE.)
          CALL PGLABEL('scan','Averaged No. of counts',INFILE)
          IF(CERR.EQ.'y')THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBIN(NSCAN,Y,ESY,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END IF
        END DO
        CALL SHOWHLP('averaged spatial profile')
C------------------------------------------------------------------------------
C Inicializamos variables
        DO I=1,NSCAN
          IFSCAN(I)=.FALSE.
        END DO
C
        DO J=1,NCHAN
          S(J)=0.
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=0.
          END DO
        END IF
C------------------------------------------------------------------------------
C Elegimos region a sumar
        WRITE(*,*)
        WRITE(*,100)'Select regions by [m]ouse of [k]eyboard (m/k) '
        CMOUSE(1:1)=READC('m','mk')
20      IF(CMOUSE.EQ.'k')THEN
21        WRITE(*,100)'Scan region (0,0=EXIT) '
          CALL READ2I('0,0',N1,N2)
          IF((N1.EQ.0).AND.(N2.EQ.0))GOTO 22
          IF((N1.LT.0).OR.(N2.GT.NSCAN).OR.(N1.GT.N2))THEN
            WRITE(*,101)'ERROR: invalid entry. Try again.'
            GOTO 21
          END IF
        ELSE
          WRITE(*,100)'Select 1st limit (q=EXIT)...'
          CALL PGBAND(0,0,0.,0.,XC,YC,CH)
          IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.(CH.EQ.'X'))THEN
            WRITE(*,*)
            GOTO 22
          END IF
          N1=NINT(XC)
          IF(N1.LT.1) N1=1
          IF(N1.GT.NSCAN) N1=NSCAN
          WRITE(*,*)N1
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSLS(2)
            CALL PGMOVE(REAL(N1),YMIN)
            CALL PGDRAW(REAL(N1),YMAX)
            IF(LCOLOR(ITERM)) CALL PGSLS(1)
          END DO
          WRITE(*,100)'Select 2nd limit (q=EXIT)...'
          CALL PGBAND(0,0,0.,0.,XC,YC,CH)
          IF((CH.EQ.'q').OR.(CH.EQ.'Q').OR.(CH.EQ.'X'))THEN
            WRITE(*,*)
            GOTO 22
          END IF
          N2=NINT(XC)
          IF(N2.LT.1) N2=1
          IF(N2.GT.NSCAN) N2=NSCAN
          WRITE(*,*)N2
          IF(N2.LT.N1)THEN
            N0=N1
            N1=N2
            N2=N1
          END IF
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGMOVE(REAL(N1),YMIN)
          CALL PGDRAW(REAL(N1),YMAX)
          CALL PGMOVE(REAL(N2),YMIN)
          CALL PGDRAW(REAL(N2),YMAX)
          CALL PGRECT(REAL(N1),REAL(N2),YMIN,YMIN+DY/100.)
          CALL PGRECT(REAL(N1),REAL(N2),YMAX,YMAX-DY/100.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
        END DO
        DO I=N1,N2
          IFSCAN(I)=.TRUE.
        END DO
        GOTO 20
C
22      N=0
        DO I=1,NSCAN
          IF(IFSCAN(I)) N=N+1
        END DO
        IF(N.EQ.0)THEN
          WRITE(*,101)'ERROR: number of scans added = 0. Try again.'
          GOTO 20
        END IF
C------------------------------------------------------------------------------
C Podemos introducir un factor multiplicativo
        WRITE(*,100)'Factor '
        FACTOR=READF('1.0')
C------------------------------------------------------------------------------
C Sumamos los nuevos espectros
        N=0
        DO I=1,NSCAN
          IF(IFSCAN(I))THEN
            N=N+1
            DO J=1,NCHAN
              S(J)=S(J)+A(J,I)*FACTOR
            END DO
            IF(CERR.EQ.'y')THEN
              DO J=1,NCHAN
                ES(J)=ES(J)+ERR(J,I)*ERR(J,I)       !el FACTOR lo metemos luego
              END DO
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
C Promediamos la suma
        DO J=1,NCHAN
          S(J)=S(J)/REAL(N)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            ES(J)=FACTOR*SQRT(ES(J))/REAL(N)
          END DO
        END IF
        WRITE(*,*)
        WRITE(*,101)'>>> File......: '//INFILE(1:TRUELEN(INFILE))
        IF(CERR.EQ.'y')THEN
          WRITE(*,101)'>>> Error file: '//ERRFILE(1:TRUELEN(ERRFILE))
        END IF
        WRITE(*,110)'>>> Number of scans added: ',N
C------------------------------------------------------------------------------
C Si es la primera imagen no dibujamos nada
        IF(LFIRSTIMAGE)THEN
          WRITE(CDUMMY,'(A2,I10)')'1,',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'First and last channel to be compared '
          CALL READ2I(CDUMMY(1:L),NC1,NC2)
C
          CALL GIVESN(NCHAN,S,ES,CERR,CMEAN1,CSNRA1)
          DO J=1,NCHAN
            FS(J)=S(J)
          END DO
          DO J=1,NCHAN
            FES(J)=ES(J)
          END DO
C
          LFIRSTIMAGE=.FALSE.
          GOTO 80
        END IF
C------------------------------------------------------------------------------
C Dibujamos
        XMIN=REAL(NC1)
        XMAX=REAL(NC2)
        DX=XMAX-XMIN
        XMIN=XMIN-DX/50.
        XMAX=XMAX+DX/50.
C
        CALL FINDMML(NCHAN,NC1,NC2,FS,YMIN,YMAX)
        IF(CERR.EQ.'y')THEN
          CALL FINDMML(NCHAN,NC1,NC2,FES,YMIN0,YMAX0)
          IF(YMIN0.LT.YMIN) YMIN=YMIN0
          IF(YMAX0.GT.YMAX) YMAX=YMAX0
        END IF
        CALL FINDMML(NCHAN,NC1,NC2,S,YMIN0,YMAX0)
        IF(YMIN0.LT.YMIN) YMIN=YMIN0
        IF(YMAX0.GT.YMAX) YMAX=YMAX0
        IF(CERR.EQ.'y')THEN
          CALL FINDMML(NCHAN,NC1,NC2,ES,YMIN0,YMAX0)
          IF(YMIN0.LT.YMIN) YMIN=YMIN0
          IF(YMAX0.GT.YMAX) YMAX=YMAX0
        END IF
C
        DY=YMAX-YMIN
        YMIN=YMIN-DY/50.
        YMAX=YMAX+DY/50.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGLABEL('channel','Averaged No. of counts',' ')
          CALL PGBIN(NCHAN,X,FS,.TRUE.)
          IF(CERR.EQ.'y')THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(3)
            CALL PGBIN(NCHAN,X,FES,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,101)'* Previous DATA'
        IF(CERR.EQ.'y')WRITE(*,101)'  [white (data), green (errors)]'
        CALL GIVESN(NCHAN,FS,FES,CERR,CMEAN1,CSNRA1)
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NCHAN,X,S,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          IF(CERR.EQ.'y')THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(4)
            CALL PGBIN(NCHAN,X,ES,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,101)'* New Spectrum'
        IF(CERR.EQ.'y')WRITE(*,101)'  [red (data), blue (errors)]'
        CALL GIVESN(NCHAN,S,ES,CERR,CMEAN2,CSNRA2)
C
        WRITE(*,*)
        CALL SHOWHLP('plot both')
        WRITE(*,100)'Press <RETURN>...'
        READ(*,*)
C------------------------------------------------------------------------------
C Comparamos Previous DATA con New Spectrum
        DO J=NC1,NC2
          SPD(J)=FS(J)-S(J)
        END DO
C Mean y Desviacion estandard
        MEAN=0.D0
        DO J=NC1,NC2
          MEAN=MEAN+SPD(J)
        END DO
        MEAN=MEAN/DBLE(NC2-NC1+1)
        SIGMA=0.D0
        DO J=NC1,NC2
          SIGMA=SIGMA+DBLE(MEAN-SPD(J))*DBLE(MEAN-SPD(J))
        END DO
        SIGMA=DSQRT(SIGMA/DBLE(NC2-NC1+1-1))
        WRITE(*,*)
        WRITE(*,101)'* Previous DATA - New Spectrum'
        WRITE(*,100)'>>> Mean value........: '
        WRITE(*,*) REAL(MEAN)
        WRITE(*,100)'>>> Standard deviation: '
        WRITE(*,*) REAL(SIGMA)
C
        CALL FINDMML(NCHAN,NC1,NC2,SPD,YMIN0,YMAX0)
        DY=YMAX0-YMIN0
        YMIN0=YMIN0-DY/50.
        YMAX0=YMAX0+DY/50.
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN0,YMAX0,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NCHAN,X,SPD,.TRUE.)
          CALL PGLABEL('channel','Averaged No. of counts',
     +     'Previous DATA - New Spectrum')
        END DO
C Ajustamos cte, recta y parabola
        DO J=NC1,NC2
          XFIT(J-NC1+1)=X(J)
          YFIT(J-NC1+1)=SPD(J)
        END DO
        NFIT=NC2-NC1+1
C
        DO NTERMS=1,3
          CALL POLFIT(XFIT,YFIT,YFIT,NFIT,NTERMS,0,COEFF,CHISQR)
          WRITE(*,*)
          DO I=1,NTERMS
            WRITE(*,'(A2,I1,A3,$)')'a(',I-1,'): '
            WRITE(*,*)COEFF(I)
          END DO
          DO J=NC1,NC2
            SPD(J)=COEFF(NTERMS)
            DO L=NTERMS-1,1,-1
              SPD(J)=SPD(J)*X(J)+COEFF(L)
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(LCOLOR(ITERM)) CALL PGSCI(NTERMS+1)
            CALL PGLINE(NCHAN,X,SPD)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END DO
        END DO
        CALL SHOWHLP('compare spectra')
C------------------------------------------------------------------------------
C Decidimos si sumamos o no el nuevo espectro
        WRITE(*,*)
        WRITE(*,100)'Add New Spectrum to Previous DATA (y/n) '
        CADD(1:1)=READC('y','yn')
        IF(CADD.EQ.'n') GOTO 80
        CALL SHOWHLP('addition result')
C------------------------------------------------------------------------------
C Hacemos una copia de Previous DATA por si queremos volver atras
        DO J=1,NCHAN
          FS0(J)=FS(J)
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            FES0(J)=FES(J)
          END DO
        END IF
C------------------------------------------------------------------------------
C Sumamos el nuevo espectro al ya acumulado
        DO J=1,NCHAN
          FS(J)=(FS(J)+S(J))/2.
        END DO
        IF(CERR.EQ.'y')THEN
          DO J=1,NCHAN
            FES(J)=SQRT(FES(J)*FES(J)+ES(J)*ES(J))/2.
          END DO
        END IF
C------------------------------------------------------------------------------
C Dibujamos el nuevo espectro
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          CALL PGIDEN_RED
          CALL PGBIN(NCHAN,X,FS,.TRUE.)
          IF(CERR.EQ.'y')THEN
            IF(LCOLOR(ITERM)) CALL PGSCI(2)
            CALL PGBIN(NCHAN,X,FES,.TRUE.)
            IF(LCOLOR(ITERM)) CALL PGSCI(1)
          END IF
          CALL PGLABEL('channel','Averaged No. of counts',
     +     'New Final DATA')
        END DO
        WRITE(*,*)
        WRITE(*,101)'* New Final DATA'
        CALL GIVESN(NCHAN,FS,FES,CERR,CMEAN3,CSNRA3)
C
        WRITE(*,*)
        WRITE(*,101)'----------------------------------------------'//
     +   '--------------------'
        WRITE(*,101)'                              (1)            '//
     +   '(2)            (3)'
        WRITE(*,100)'>>> Mean spectrum value: '
        WRITE(*,101)CMEAN1//CMEAN2//CMEAN3
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'>>> <S/N ratio>/pixel..: '
          WRITE(*,101)CSNRA1//CSNRA2//CSNRA3
        END IF
        WRITE(*,101)'----------------------------------------------'//
     +   '--------------------'
        WRITE(*,101)'(1): Previous DATA'
        WRITE(*,101)'(2): New Spectrum'//' ['//
     +   INFILE(1:TRUELEN(INFILE))//']'
        WRITE(*,101)'(3): New Final DATA [(1)+(2)]'
        WRITE(*,101)'----------------------------------------------'//
     +   '--------------------'
C------------------------------------------------------------------------------
C Pedimos confirmacion del nuevo resultado
        CALL SHOWHLP('final decision')
        WRITE(*,100)'Accept New Final Spectrum (y/n) '
        CCONF(1:1)=READC('y','yn')
        IF(CCONF.EQ.'n')THEN
          WRITE(*,100)'Restoring Previous DATA...'
          DO J=1,NCHAN
            FS(J)=FS0(J)
          END DO
          IF(CERR.EQ.'y')THEN
            DO J=1,NCHAN
              FES(J)=FES0(J)
            END DO
          END IF
          WRITE(*,101)'  ...OK!'
        END IF
C------------------------------------------------------------------------------
80      WRITE(*,*)
        CALL SHOWHLP('more files')
        WRITE(*,100)'Are you using a new file (y/n) '
        CNEW(1:1)=READC('y','yn')
        IF(CNEW.EQ.'y') GOTO 10
C------------------------------------------------------------------------------
C Cerramos la salida grafica
        CALL PGEND
C------------------------------------------------------------------------------
C Preguntamos si salvamos fichero o no
        WRITE(*,*)
        WRITE(*,100)'Save New Final DATA into file (y/n) '
        CSAVE(1:1)=READC('y','yn')
        IF(CSAVE.EQ.'n') STOP
C------------------------------------------------------------------------------
C Salvamos fichero
        WRITE(*,*)
        CALL SHOWHLP('expansion')
        WRITE(*,100)'Expand final spectrum (y/n) '
        CEXPAND(1:1)=READC('n','yn')
        IF(CEXPAND.EQ.'y')THEN
          WRITE(*,100)'NSCAN in output file'
          NSCAN=READI('@')
        ELSE
          NSCAN=1
        END IF
        WRITE(*,100)'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (FS(J),J=1,NCHAN)
        END DO
        CLOSE(30)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(31,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(31) (FES(J),J=1,NCHAN)
          END DO
          CLOSE(31)
        END IF
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Proporciona la media y (si CERR='y') la relacion senhal ruido
C Las variables de texto CMEAN y CSNRA devuelven (en modo caracter), el
C valor medio y (si procede) la relacion senhal/ruido
        SUBROUTINE GIVESN(N,S,ES,CERR,CMEAN,CSNRA)
        IMPLICIT NONE
C
        INTEGER N
        REAL S(N),ES(N)
        CHARACTER*1 CERR
        CHARACTER*(*) CMEAN
        CHARACTER*(*) CSNRA
C
        INTEGER NC1,NC2
        INTEGER J,NSN
        DOUBLE PRECISION MEAN,SN
C
        COMMON/BLKNC/NC1,NC2
C------------------------------------------------------------------------------
        WRITE(CMEAN,*) 0.0
        WRITE(CSNRA,*) 0.0
C
        MEAN=0.D0
        DO J=NC1,NC2
          MEAN=MEAN+DBLE(S(J))
        END DO
        MEAN=MEAN/DBLE(NC2-NC1+1)
        WRITE(*,100)'>>> Mean spectrum value: '
        WRITE(*,*) REAL(MEAN)
        WRITE(CMEAN,*) REAL(MEAN)
C
        IF(CERR.EQ.'y')THEN
          NSN=0
          SN=0.D0
          DO J=NC1,NC2
            IF((S(J).GT.0.).AND.(ES(J).GT.0.))THEN
              NSN=NSN+1
              SN=SN+DBLE(S(J)/ES(J))
            END IF
          END DO
          IF(NSN.GT.0)THEN
            SN=SN/DBLE(NSN)
            WRITE(*,100)'>>> <S/N ratio>/pixel..: '
            WRITE(*,*)REAL(SN)
            WRITE(CSNRA,*) REAL(SN)
          END IF
          IF(NSN.NE.(NC2-NC1+1))THEN
            WRITE(*,101)'>>> WARNING: data and/or errors .le. 0 '//
     +       'have been found.'
          END IF
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
