C------------------------------------------------------------------------------
C Version 13-April-1999                                           file: imath.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: imath
C Classification: arithmetic & manipulations
C Description: Performs arithmetic with images, spectra, spatial cross sections
C and constants.
C
Comment
C
C Realiza sobre una imagen operaciones algebraicas elementales (+-/*) con
C otras imagenes, espectros, cortes espaciales y constantes. Tambien puede
C calcular la imagen de errores correspondiente.
C
        PROGRAM IMATH
        IMPLICIT NONE
C Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J,L
        INTEGER ITS,IOPER
        INTEGER NSCAN2,NCHAN2
        REAL A(NCMAX,NSMAX),ERR(NCMAX,NSMAX),S(NCMAX),ES(NCMAX)
        REAL STWV2,DISP2
        REAL CTE,ECTE,CTE2,ECTE2
        REAL FNULL,ERRFNULL
        CHARACTER*1 CERR
        CHARACTER*50 CDUMMY
        CHARACTER*75 INFILE,ERRFILE,OUTFILE
        LOGICAL LAVOID
C
C Si se cambia redlib.inc, tambien hay que introducir los cambios aqui
        REAL AIRMASS2,TIMEXPOS2
        CHARACTER*255 OBJECT2,FITSFILE2,COMMENT2
C------------------------------------------------------------------------------
        THISPROGRAM='imath'
        CALL WELCOME('13-April-1999')
C avoid compilation warnings
        FNULL=0.0
        ERRFNULL=0.0
        CTE=0.0
        CTE2=0.0
        ECTE2=0.0
C------------------------------------------------------------------------------
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        AIRMASS2=AIRMASS      !conservamos las variables para el fichero output
        TIMEXPOS2=TIMEXPOS
        OBJECT2=OBJECT
        FITSFILE2=FITSFILE
        COMMENT2=COMMENT
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Input error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(21,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(21) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(21)
        END IF
C
        WRITE(*,*)
        WRITE(*,101)'OPERATION:'
        WRITE(*,101)'  1 = add file/cte          2 = subtract file/cte'
        WRITE(*,101)'  3 = multiply by file/cte  4 = divide by file/cte'
        WRITE(*,100)'  5 = raise to a power'
        WRITE(*,101)' (be careful with negative or null data)'
        WRITE(*,101)'  6 = divide by file/cte avoiding division by 0'
ccc10      WRITE(*,100)'Option '
        WRITE(*,100)'Option '
        IOPER=READILIM('@',1,6)
        IF(IOPER.EQ.6)THEN
          WRITE(*,100) 'Real number to be obtained in case of '
          WRITE(*,100) 'division by zero '
          FNULL=READF('0.0')
          IF(CERR.EQ.'y')THEN
            WRITE(*,100) 'Error for this number '
            ERRFNULL=READF('0.0')
          END IF
          IOPER=4
          LAVOID=.TRUE.
        ELSE
          LAVOID=.FALSE.
        END IF
C
        WRITE(*,*)
ccc12      WRITE(*,101)'2ND FILE WILL BE A:'
        WRITE(*,101)'2ND FILE WILL BE A:'
        WRITE(*,100)'  0 = frame        --> '
        WRITE(CDUMMY,*)NSCAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'NSCAN: '//CDUMMY(1:L)
        WRITE(CDUMMY,*)NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)', NCHAN: '//CDUMMY(1:L)
        WRITE(*,100)'  1 = spectrum     --> '
        WRITE(CDUMMY,*)NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'NCHAN: '//CDUMMY(1:L)
        WRITE(*,100)'  2 = spatial form --> '
        WRITE(CDUMMY,*)NSCAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,101)'NSCAN: '//CDUMMY(1:L)
        WRITE(*,101)'  3 = constant     --> (by keyboard)'
        WRITE(*,100)'Option '
        ITS=READILIM('@',0,3)
C
        WRITE(*,*)
        IF(ITS.EQ.0)THEN
          WRITE(*,101)'Enter 2nd image:'
        ELSEIF(ITS.EQ.1)THEN
          WRITE(*,101)'Enter spectrum:'
        ELSEIF(ITS.EQ.2)THEN
          WRITE(*,101)'Enter spatial form: '
        ELSE
          WRITE(*,100)'Constant'
          CTE=READF('@')
        END IF
        IF(ITS.NE.3)THEN
          WRITE(*,100)'Input file name'
          INFILE=INFILEX(30,'@',NSCAN2,NCHAN2,STWV2,DISP2,1,.FALSE.)
        END IF
C
        IF(CERR.EQ.'y')THEN
          IF(ITS.EQ.0)THEN
            WRITE(*,101)'Enter 2nd image error:'
          ELSEIF(ITS.EQ.1)THEN
            WRITE(*,101)'Enter spectrum error: '
          ELSEIF(ITS.EQ.2)THEN
            WRITE(*,101)'Enter spatial form error: '
          ELSE
            WRITE(*,100)'Error in constant'
            ECTE=READF('@')
            ECTE2=ECTE*ECTE
            CTE2=CTE*CTE
          END IF
          IF(ITS.NE.3)THEN
            WRITE(*,100)'Input error file name '
            CALL GUESSEF(INFILE,ERRFILE)
            ERRFILE=INFILEX(31,ERRFILE,NSCAN2,NCHAN2,STWV2,DISP2,21,
     +       .TRUE.)                                                  !match
          END IF
        END IF
C
        IF(ITS.NE.3)THEN
          IF((STWV.NE.STWV2).OR.(DISP.NE.DISP2))THEN
            WRITE(*,101)'WARNING: STWV and/or DISP of 1st and 2nd '//
     +       'files are different.'
            WRITE(*,100)'Press <RETURN>...'
            READ(*,*)
          END IF
        END IF
C------------------------------------------------------------------------------
        IF(ITS.EQ.0)THEN                                     !2nd FRAME = IMAGE
          IF((NSCAN.NE.NSCAN2).OR.(NCHAN.NE.NCHAN2))THEN
            WRITE(*,101)'FATAL ERROR: 1st and 2nd files do not have '//
     +       'the same size.'
            CLOSE(30)
            IF(CERR.EQ.'y') CLOSE(31)
            STOP
          END IF
          DO I=1,NSCAN
            READ(30) (S(J),J=1,NCHAN)
            IF(CERR.EQ.'y') READ(31) (ES(J),J=1,NCHAN)
            IF(IOPER.EQ.1)THEN                                  !suma
              IF(CERR.EQ.'y')THEN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ES(J)*ES(J))
                END DO
              END IF
              DO J=1,NCHAN
                A(J,I)=A(J,I)+S(J)
              END DO
            ELSEIF(IOPER.EQ.2)THEN                              !resta
              IF(CERR.EQ.'y')THEN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ES(J)*ES(J))
                END DO
              END IF
              DO J=1,NCHAN
                A(J,I)=A(J,I)-S(J)
              END DO
            ELSEIF(IOPER.EQ.3)THEN                              !multiplicacion
              IF(CERR.EQ.'y')THEN
                DO J=1,NCHAN
                  ERR(J,I)=A(J,I)*A(J,I)*ES(J)*ES(J)+
     +             S(J)*S(J)*ERR(J,I)*ERR(J,I)
                  ERR(J,I)=SQRT(ERR(J,I))
                END DO
              END IF
              DO J=1,NCHAN
                A(J,I)=A(J,I)*S(J)
              END DO
            ELSEIF(IOPER.EQ.4)THEN                              !division
              IF(CERR.EQ.'y')THEN
                IF(LAVOID)THEN
                  DO J=1,NCHAN
                    IF(S(J).EQ.0.0)THEN
                      ERR(J,I)=ERRFNULL
                    ELSE
                      ERR(J,I)=A(J,I)*A(J,I)*ES(J)*ES(J)+
     +                 S(J)*S(J)*ERR(J,I)*ERR(J,I)
                      ERR(J,I)=SQRT(ERR(J,I))/(S(J)*S(J))
                    END IF
                  END DO
                ELSE
                  DO J=1,NCHAN
                    ERR(J,I)=A(J,I)*A(J,I)*ES(J)*ES(J)+
     +               S(J)*S(J)*ERR(J,I)*ERR(J,I)
                    ERR(J,I)=SQRT(ERR(J,I))/(S(J)*S(J))
                  END DO
                END IF
              END IF
              IF(LAVOID)THEN
                DO J=1,NCHAN
                  IF(S(J).EQ.0.0)THEN
                    A(J,I)=FNULL
                  ELSE
                    A(J,I)=A(J,I)/S(J)
                  END IF
                END DO
              ELSE
                DO J=1,NCHAN
                  A(J,I)=A(J,I)/S(J)
                END DO
              END IF
            ELSEIF(IOPER.EQ.5)THEN                  !elevamos datos a potencias
              IF(CERR.EQ.'y')THEN
                DO J=1,NCHAN
                  ERR(J,I)=S(J)*S(J)*A(J,I)**(2.*S(J)-2.)*
     +             ERR(J,I)*ERR(J,I)+
     +             A(J,I)**(2.*S(J))*ALOG(A(J,I))*ALOG(A(J,I))*
     +             ES(J)*ES(J)
                  ERR(J,I)=SQRT(ERR(J,I))
                  END DO
              END IF
              DO J=1,NCHAN
                A(J,I)=A(J,I)**S(J)
              END DO
            END IF
          END DO
          CLOSE(30)
          IF(CERR.EQ.'y') CLOSE(31)
C------------------------------------------------------------------------------
        ELSEIF(ITS.EQ.1)THEN                              !2nd FRAME = SPECTRUM
          IF(NCHAN.NE.NCHAN2)THEN
            WRITE(*,101)'ERROR: spectrum size does not match image'//
     +       ' size.'
            CLOSE(30)
            IF(CERR.EQ.'y') CLOSE(31)
            STOP
          END IF
          IF(NSCAN2.NE.1)THEN
            WRITE(*,100)'WARNING: NSCAN.GT.1 in last frame'
            WRITE(*,101)' (only 1st scan will be read and employed)'
          END IF
          READ(30) (S(J),J=1,NCHAN)
          CLOSE(30)
          IF(CERR.EQ.'y')THEN
            READ(31) (ES(J),J=1,NCHAN)
            CLOSE(31)
          END IF
          IF(IOPER.EQ.1)THEN                                    !suma
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ES(J)*ES(J))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)+S(J)
              END DO
            END DO
          ELSEIF(IOPER.EQ.2)THEN                                !resta
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ES(J)*ES(J))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)-S(J)
              END DO
            END DO
          ELSEIF(IOPER.EQ.3)THEN                                !multiplicacion
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=A(J,I)*A(J,I)*ES(J)*ES(J)+
     +             S(J)*S(J)*ERR(J,I)*ERR(J,I)
                  ERR(J,I)=SQRT(ERR(J,I))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)*S(J)
              END DO
            END DO
          ELSEIF(IOPER.EQ.4)THEN                                !division
            IF(CERR.EQ.'y')THEN
              IF(LAVOID)THEN
                DO I=1,NSCAN
                  DO J=1,NCHAN
                    IF(S(J).EQ.0.0)THEN
                      ERR(J,I)=ERRFNULL
                    ELSE
                      ERR(J,I)=A(J,I)*A(J,I)*ES(J)*ES(J)+
     +                 S(J)*S(J)*ERR(J,I)*ERR(J,I)
                      ERR(J,I)=SQRT(ERR(J,I))/(S(J)*S(J))
                    END IF
                  END DO
                END DO
              ELSE
                DO I=1,NSCAN
                  DO J=1,NCHAN
                    ERR(J,I)=A(J,I)*A(J,I)*ES(J)*ES(J)+
     +               S(J)*S(J)*ERR(J,I)*ERR(J,I)
                    ERR(J,I)=SQRT(ERR(J,I))/(S(J)*S(J))
                  END DO
                END DO
              END IF
            END IF
            IF(LAVOID)THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  IF(S(J).EQ.0.0)THEN
                    A(J,I)=FNULL
                  ELSE
                    A(J,I)=A(J,I)/S(J)
                  END IF
                END DO
              END DO
            ELSE
              DO I=1,NSCAN
                DO J=1,NCHAN
                  A(J,I)=A(J,I)/S(J)
                END DO
              END DO
            END IF
          ELSEIF(IOPER.EQ.5)THEN                    !elevamos datos a potencias
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=S(J)*S(J)*A(J,I)**(2.*S(J)-2.)*
     +             ERR(J,I)*ERR(J,I)+
     +             A(J,I)**(2.*S(J))*ALOG(A(J,I))*ALOG(A(J,I))*
     +             ES(J)*ES(J)
                  ERR(J,I)=SQRT(ERR(J,I))
                  END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)**S(J)
              END DO
            END DO
          END IF
C------------------------------------------------------------------------------
        ELSEIF(ITS.EQ.2)THEN                          !2nd FRAME = SPATIAL FORM
          IF(NSCAN.NE.NSCAN2)THEN
            WRITE(*,101)'ERROR: spatial form size does not math'//
     +       ' image size.'
            CLOSE(30)
            IF(CERR.EQ.'y') CLOSE(31)
            STOP
          END IF
          IF(NCHAN2.NE.1)THEN
            WRITE(*,100)'WARNING: NCHAN.GT.1 in last frame'
            WRITE(*,101)' (only 1st channel will be read and employed)'
          END IF
          DO I=1,NSCAN
            READ(30) S(I)
          END DO
          CLOSE(30)
          IF(CERR.EQ.'y')THEN
            DO I=1,NSCAN
              READ(31) ES(I)
            END DO
            CLOSE(31)
          END IF
          IF(IOPER.EQ.1)THEN                                    !suma
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ES(I)*ES(I))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)+S(I)
              END DO
            END DO
          ELSEIF(IOPER.EQ.2)THEN                                !resta
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ES(I)*ES(I))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)-S(I)
              END DO
            END DO
          ELSEIF(IOPER.EQ.3)THEN                                !multiplicacion
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=A(J,I)*A(J,I)*ES(I)*ES(I)+
     +             S(I)*S(I)*ERR(J,I)*ERR(J,I)
                  ERR(J,I)=SQRT(ERR(J,I))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)*S(I)
              END DO
            END DO
          ELSEIF(IOPER.EQ.4)THEN                                !division
            IF(CERR.EQ.'y')THEN
              IF(LAVOID)THEN
                DO I=1,NSCAN
                  DO J=1,NCHAN
                    IF(S(I).EQ.0.0)THEN
                      ERR(J,I)=ERRFNULL
                    ELSE
                      ERR(J,I)=A(J,I)*A(J,I)*ES(I)*ES(I)+
     +                 S(I)*S(I)*ERR(J,I)*ERR(J,I)
                      ERR(J,I)=SQRT(ERR(J,I))/(S(I)*S(I))
                    END IF
                  END DO
                END DO
              ELSE
                DO I=1,NSCAN
                  DO J=1,NCHAN
                    ERR(J,I)=A(J,I)*A(J,I)*ES(I)*ES(I)+
     +               S(I)*S(I)*ERR(J,I)*ERR(J,I)
                    ERR(J,I)=SQRT(ERR(J,I))/(S(I)*S(I))
                  END DO
                END DO
              END IF
            END IF
            IF(LAVOID)THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  IF(S(I).EQ.0.0)THEN
                    A(J,I)=FNULL
                  ELSE
                    A(J,I)=A(J,I)/S(I)
                  END IF
                END DO
              END DO
            ELSE
              DO I=1,NSCAN
                DO J=1,NCHAN
                  A(J,I)=A(J,I)/S(I)
                END DO
              END DO
            END IF
          ELSEIF(IOPER.EQ.5)THEN                    !elevamos datos a potencias
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=S(I)*S(I)*A(J,I)**(2.*S(I)-2.)*
     +             ERR(J,I)*ERR(J,I)+
     +             A(J,I)**(2.*S(I))*ALOG(A(J,I))*ALOG(A(J,I))*
     +             ES(I)*ES(I)
                  ERR(J,I)=SQRT(ERR(J,I))
                  END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)**S(I)
              END DO
            END DO
          END IF
C------------------------------------------------------------------------------
        ELSE                                                   !2nd FRAME = Cte
          IF(IOPER.EQ.1)THEN                                    !suma
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ECTE2)
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)+CTE
              END DO
            END DO
          ELSEIF(IOPER.EQ.2)THEN                                !resta
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=SQRT(ERR(J,I)*ERR(J,I)+ECTE2)
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)-CTE
              END DO
            END DO
          ELSEIF(IOPER.EQ.3)THEN                                !multiplicacion
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=A(J,I)*A(J,I)*ECTE2+CTE2*ERR(J,I)*ERR(J,I)
                  ERR(J,I)=SQRT(ERR(J,I))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)*CTE
              END DO
            END DO
          ELSEIF(IOPER.EQ.4)THEN                                !division
            IF(CERR.EQ.'y')THEN
              IF(LAVOID)THEN
                DO I=1,NSCAN
                  DO J=1,NCHAN
                    IF(CTE2.EQ.0.0)THEN
                      ERR(J,I)=ERRFNULL
                    ELSE
                      ERR(J,I)=A(J,I)*A(J,I)*ECTE2+
     +                         CTE2*ERR(J,I)*ERR(J,I)
                      ERR(J,I)=SQRT(ERR(J,I))/(CTE2)
                    END IF
                  END DO
                END DO
              ELSE
                DO I=1,NSCAN
                  DO J=1,NCHAN
                    ERR(J,I)=A(J,I)*A(J,I)*ECTE2+CTE2*ERR(J,I)*ERR(J,I)
                    ERR(J,I)=SQRT(ERR(J,I))/(CTE2)
                  END DO
                END DO
              END IF
            END IF
            IF(LAVOID)THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  IF(CTE.EQ.0.0)THEN
                    A(J,I)=FNULL
                  ELSE
                    A(J,I)=A(J,I)/CTE
                  END IF
                END DO
              END DO
            ELSE
              DO I=1,NSCAN
                DO J=1,NCHAN
                  A(J,I)=A(J,I)/CTE
                END DO
              END DO
            END IF
          ELSEIF(IOPER.EQ.5)THEN                    !elevamos datos a potencias
            IF(CERR.EQ.'y')THEN
              DO I=1,NSCAN
                DO J=1,NCHAN
                  ERR(J,I)=CTE2*A(J,I)**(2.*CTE-2.)*ERR(J,I)*ERR(J,I)+
     +             A(J,I)**(2.*CTE)*ALOG(A(J,I))*ALOG(A(J,I))*ECTE2
                  ERR(J,I)=SQRT(ERR(J,I))
                END DO
              END DO
            END IF
            DO I=1,NSCAN
              DO J=1,NCHAN
                A(J,I)=A(J,I)**CTE
              END DO
            END DO
          END IF
        END IF
C------------------------------------------------------------------------------
C salvamos fichero final
        WRITE(*,*)
        WRITE(*,100)'Output file name'
        AIRMASS=AIRMASS2     !conservamos las variables para el fichero output
        TIMEXPOS=TIMEXPOS2
        OBJECT=OBJECT2
        FITSFILE=FITSFILE2
        COMMENT=COMMENT2
        OUTFILE=OUTFILEX(50,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(50) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(50)
C
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          OUTFILE=OUTFILEX(51,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(51) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(51)
        END IF
C
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
        END
