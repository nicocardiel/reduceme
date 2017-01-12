C------------------------------------------------------------------------------
C Version 22-October-2001                                      file: cleanest.f
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
C Program: cleanest
C Classification: cosmic rays
C Description: Automatic and interactive removal of cosmic rays.
C
Comment
C
C Eliminacion manual y/o automatica de rayos cosmicos
C------------------------------------------------------------------------------
C Version 7-Jun-1994
C Version 9-Jun-1994
C Revisited 24-Jun-1994: Option 'r' removes all the crosses
C Version 4-Jan-1995: Unix Version. Buttons added.
C Version 10-Octubre-1995: 3D-plots
C Version 13-Octubre-1995: Rotate 3D-plots
C Version 17-Abril-1996: Incorpora plotsp3d
C Version 8-Mayo-1996: Uso de imagenes auxiliares en paralelo
C Version 20-Junio-1996: eliminacion de la forma x/y promedio en el cuadrado
C  de busqueda para refinar la busqueda de C.R.'s
C Version 1-Octubre-1996: ajusta niveles medios de Aux Frames
C Version 21-Noviembre-1996: interpolamos imagen de errores simultaneamente
C Version 29-Octubre-1997: Opcion Top100 incorporada
C Version 22-Octubre-2001: corregimos error en subrutina PLOT3DBARS
C------------------------------------------------------------------------------
        PROGRAM CLEANEST
C------------------------------------------------------------------------------
C
C Este programa permite eliminar rayos cosmicos combinando busqueda
C automatica con busqueda manual. El usuario puede definir facilmente
C que regiones deben ser limpiadas y cuales no, etc.
C
C A continuacion se describen brevemente la estructura y opciones 
C principales del programa.
C
C El programa esta preparado para leer imagenes de tamanho variable
C dentro de los limites definidos por los parametros NSMAX (numero maximo
C de scans) y NCMAX (numero maximo de canales).
C
C El usuario define inicialmente el sistema de busqueda que el programa
C utilizara para detectar rayos cosmicos. Hay tres sistemas posibles (se han
C definidos los mismos tres sistemas que utiliza el comando BCLEAN de FIGARO).
C Estos sistemas vienen definidos por:
C SIGCR:  numero de veces SIGMA que debe exceder un pixel el valor medio MEAN
C         calculado por el programa en una ventana de pixels centrada en el
C         pixel en examen.
C FACTCR: numero de veces que un pixel debe exceder el valor medio MEAN 
C         para ser considerado rayo cosmico.
C MINVCR: valor minimo que debe exceder un pixel el valor medio MEAN para 
C         ser considerado rayo cosmico.
C
C Tambien hay que definir el tamanho del SEARCHING SQUARE como el tamanho 
C del cuadrado (centrado en el pixel en cuestion --salvo en los bordes, donde 
C el cuadrado se desplaza para mantener constante el numero de pixels--) 
C donde el programa mide, inicialmente, MEAN y SIGMA, para la busqueda 
C inicial de rayos cosmicos.
C
C Los tres metodos son compatibles simultaneamente. Es decir, si el usuario
C define los tres factores, el programa detectara como rayo cosmico 
C cualquier pixel que verifique uno cualquiera de los tres metodos de busqueda.
C
C Existe una cuestion fundamental que puede producir un poco de confunsion
C al usuario y que merece la pena ser tratada aqui con un poco de detalle.
C Una cosa es el sistema de busqueda de rayo cosmico y otra es el sistema 
C que el programa utiliza, una vez detectado el rayo cosmico, para determinar 
C si los pixels de alrededor son tambien (o no) rayo cosmico. El hecho de 
C definir un sistema de busqueda inicial diferente al sistema de busqueda 
C final (que vamos a llamar sistema de busqueda fino) se debe a que una vez 
C detectado un rayo cosmico resulta mas efectivo variar las condiciones de 
C busqueda para que el programa sea mas "inteligente" y pueda discernir con 
C mas facilidad que pixels forman parte del posible rayo cosmico.
C
C Por ello, una vez detectado el rayo cosmico, el programa se olvida de las 
C condiciones de busqueda iniciales (cualquiera que haya sido el metodo(s) 
C elegido), para pasar a una busqueda en la cual ya solo se tiene en cuenta 
C el valor de SIGMA y MEAN, aunque en este caso el tamanho del cuadrado de 
C busqueda ya no tiene porque coincidir con el tamanho del cuadrado donde 
C realmente se mide el valor medio y la sigma.
C
C En el sistema de busqueda fino el tamanho de la ventana sobre la cual se 
C miden MEAN y SIGMA, viene dado por NPP, el tamanho de la ventana que 
C realmente se dibuja cuando el programa busca rayos cosmicos (que puede ser 
C diferente al tamnanho del SEARCHING SQUARE). Para resumir un poco la 
C situacion, y tratar de no enfollonar mas la explicacion, vamos a representar 
C en un esquema como es el procedimiento de busqueda desde el principio:
C 
C 1) El usuario define el sistema de busqueda (SIGCR, FACTCR, MINVCR) y el 
C tamanho del SEARCHING SQUARE.
C 
C 2) Segun la eleccion anterior, el programa detectara como rayo cosmico 
C aquellos pixels que verifiquen alguna de las condiciones de busqueda, 
C teniendo en cuenta que MEAN y SIGMA se calculan con todos los pixels del 
C SEARCHING SQUARE que no estan en el mismo scan y canal que el pixel en 
C estudio (eliminando los pixels del mismo scan y canal que el pixel 
C estudiado podemos eliminar en los calculos posibles pixels tambien 
C afectados por rayo cosmico).
C 
C 3) Si se encuentra el rayo cosmico, el programa hace un primer plot con los 
C cortes en +3 SIGMA y -3 SIGMA alrededor de MEAN. Este plot representa un 
C cuadrado de NPP x NPP pixels, (aunque insisto, el tamanho del SEARCHING 
C SQUARE sobre el cual se han medido MEAN y SIGMA haya sido mayor o 
C menor.)
C 
C 4) Ahora, usando un cuadrado de tamanho NPP x NPP, el programa ordenada 
C todos los pixels de dicho cuadrado segun su senhal. Utilizando entonces una 
C fraccion de estos pixels con la senhal comprendida entre PERCENT1 y PERCENT2
C (ver opcion # Special Options), se calcula un nuevo valor de MEAN y SIGMA. 
C El programa suele calcular un valor medio mucho mas ajustado al 
C nivel de fondo, con lo cual se han ignorado todos los pixels del rayo 
C cosmico (en muchos casos incluso la senhal de lineas de cielo), y posibles
C pixels "muertos".
C
C 5) Un nuevo plot NPP x NPP aparece ahora con los cortes en +SIGTHRESHOLD y
C -SIGTHRESHOLD veces SIGMA alrededor de MEAN.
C 
C 6) El programa presenhaliza ahora los pixels que superen MEAN en 
C SIGTHRESHOLD veces SIGMA. Normalmente el valor de SIGTHRESHOLD es un 50% el 
C valor de SIGCR (en caso de haber utilizado este metodo). De esta forma el 
C programa puede detectar no solo los pixels que forman parte mas claramente de
C un rayo cosmico, sino ademas aquellos afectados en menor medida, 
C pero que tambien tienen que interpolarse.
C 
C 7) El usuario puede entonces desmarcar alguno de los pixels presenhalizados
C por el programa o/y anhadir alguno mas.
C 
C 8) Entonces basta con interpolar diciendole al programa si la interpolacion 
C se realiza siguiendo el eje X o el eje Y. La opcion N (number of pixels to 
C be employed at both sides) determina el numero de pixels a cada lado que el 
C programa utilizara para determinar el valor a sustituir. Este se calcula
C mediante el ajuste a un polinomio de grado NDP-1 (NDP=numero de terminos
C del polinomio). Este polinomio se calcula ajustando los NFIT puntos a cada
C lado del conjunto de pixels marcados (consecutivamente). Hay que indicar
C que no es necesario marcar todos los pixels en horizontal o vertical para
C que estos se consideren consecutivos. El programa determina el primero y
C el ultimo (en vertical u horizontal) y considera que todos los intermedios
C estan marcados (aunque no sea asi de hecho). Esto facilita la eliminacion
C de rayos cosmicos grandes, dado que no hay que marcar todos los pixels
C afectados.
C 
C El programa permite asimismo ver la imagen completa, hacer Zoom y buscar 
C con el cursor rayos cosmicos siguiendo la filosolia del comando CLEAN de 
C FIGARO.
C
C Facilidades adicionales pueden encontrarse en los diferentes menus que el
C programa ofrece.
C
C Hemos introducido asimismo la posibilidad de interpolar simultaneamente en
C las imagenes de errores. El programa simplemente repite la interpolacion
C realizada sobre la imagen de datos (usando los mismos pixels y grado del
C polinomio).
C
C------------------------------------------------------------------------------
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NPPMAX
        INTEGER NBOTONES
        INTEGER NAUXMAX                   !numero maximo de imagenes auxiliares
        INTEGER NTOPMAX
        PARAMETER(NPPMAX=101)         !si se cambia hacerlo en todos los sitios
        PARAMETER(NBOTONES=18)
        PARAMETER (NAUXMAX=3)         !si se cambia hacerlo en todos los sitios
        PARAMETER (NTOPMAX=1000)                     !numero maximo para TOP100
C
        INTEGER I,J,K,L,L1,L2
        INTEGER NPP
        INTEGER NS,NC
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NC1P,NC2P,NS1P,NS2P
        INTEGER NWIN
        INTEGER NFIT
        INTEGER NDP  !grado del polinomio+1 en los ajustes para eliminar C.R.'s
        INTEGER IMENU
        INTEGER NS1YES,NS2YES,NC1YES,NC2YES,NS0YES,NC0YES
        INTEGER CAUTOMODE
        INTEGER NTOTCR
        INTEGER BMODE(NBOTONES),BMODE2(NBOTONES)
        INTEGER NB
        INTEGER NCRSAVED
        INTEGER NTYPEPLOT
        INTEGER NAUXFRAME
        INTEGER ANGLE3D
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NTOP,ITOP(NTOPMAX),JTOP(NTOPMAX)
        REAL A(NCMAX,NSMAX)
        REAL ERR(NCMAX,NSMAX)
        REAL AUXFRAME(NCMAX,NSMAX,NAUXMAX)
        REAL MEAN,SIGMA,DEV
        REAL SIGCR,FACTCR,MINVCR,SIGCRAUX
        REAL PERCENT1,PERCENT2
        REAL SIGTHRESHOLD
        REAL RQSKY
        REAL XC,YC
        REAL FACTOR
        REAL FIXED_FG,FIXED_BG
        REAL BG,FG,TR(6),MINVAL,MAXVAL
        REAL ATOP(NTOPMAX)
        REAL AREPLACE(NPPMAX,NPPMAX)
        CHARACTER*1 CH,CDEF,CSURE,CPINTAAUTO,CCROS,CSKIPCH,CMOUSE
        CHARACTER*1 CCLEAN2,CFIX_FGBG,CERR
        CHARACTER*20 LABEL(NBOTONES)
        CHARACTER*20 LABEL2(NBOTONES)
        CHARACTER*50 CDUMMY
        CHARACTER*75 AUXFILE(NAUXMAX),ERRFILE
        CHARACTER*80 NAMFIL,HISTFILE
        LOGICAL TOP100CLEAN(NCMAX,NSMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        LOGICAL CRFOUND
        LOGICAL LNEXT
        LOGICAL CLEANAUTO
        LOGICAL SCANYES(NSMAX),CHANYES(NCMAX)
        LOGICAL INSIDE2
        LOGICAL LHISTFILE
        LOGICAL LREMOVEX,LREMOVEY
        LOGICAL LEXIT
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLKERR/ERR
        COMMON/BLK1/MEAN,SIGMA  
        COMMON/BLK2/NS,NC,NWIN
        COMMON/BLK3A/LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        COMMON/BLK3B/SIGCR,FACTCR,MINVCR,SIGCRAUX
        COMMON/BLK5/NFIT,NPP,NDP
        COMMON/BLK6/PERCENT1,PERCENT2
        COMMON/BLK7A/SIGTHRESHOLD
        COMMON/BLK7B/RQSKY
        COMMON/BLK9/CLEANAUTO,CAUTOMODE,CCROS
        COMMON/BLK10/NTOTCR
        COMMON/BLK11/LABEL2,BMODE2
        COMMON/BLK12/NAMFIL
        COMMON/BLK13/LHISTFILE
        COMMON/BLK14/NCRSAVED
        COMMON/BLK15/NTYPEPLOT
        COMMON/BLK16/FACTOR
        COMMON/BLK17A/NAUXFRAME
        COMMON/BLK17B/AUXFRAME
        COMMON/BLK17C/AUXFILE
        COMMON/BLK17D/STWV,DISP
        COMMON/BLK18/LREMOVEX,LREMOVEY
        COMMON/BLK19/CHANYES,SCANYES
        COMMON/BLK20A/CFIX_FGBG
        COMMON/BLK20B/FIXED_FG,FIXED_BG
        COMMON/BLK21/ANGLE3D
        COMMON/BLK22/CERR
        COMMON/BLK23/AREPLACE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
        COMMON/BLK1TOP100/ITOP,JTOP
        COMMON/BLK2TOP100/ATOP
        COMMON/BLK3TOP100/TOP100CLEAN
C------------------------------------------------------------------------------
        DATA (LABEL(I),I=1,NBOTONES)/
     +   '[s]tart','[r]egion','[w]indow','[a]utomatic','[l]ook',
     +   '[o]ptions',
     +   'sa[v]e','histo->[1]','histo->[2]','[t]op100','plotsp3[d]',
     +   '[q]uit',
     +   ' ',' ',' ',' ',' ',' '/
        DATA (BMODE(I),I=1,NBOTONES) / 0, 0, 0, 0, 0, 0,
     +                                 0, 0, 0, 0, 0, 0,
     +                                -1,-1,-1,-1,-1,-1/
        DATA (LABEL2(I),I=1,NBOTONES)/
     +   '[x](interp)','[y](interp)','surf[a]ce','repla[z]e',
     +   '[s](bg/fg)','[n]=','[r](crosses)','[t](test)',
     +   '[u](undo)','[#](options)','[c](cont.)','[e](Exit)',
     +   '3D([+]90\uo\d)','3D([-]90\uo\d)','3D(auto)',
     +   'Aux.[f]rame','s[k]ip x','sk[i]p y'/
        DATA (BMODE2(I),I=1,NBOTONES)/ 0, 0, 0, 0, 0, 0,
     +                                 0, 0, 3, 0, 0, 0,
     +                                 3, 3, 3, 3, 0, 0/
C------------------------------------------------------------------------------
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C------------------------------------------------------------------------------
        THISPROGRAM='cleanest'
        CALL WELCOME('29-November-1997')
        CALL SHOWHLP('explanation')
C------------------------------------------------------------------------------
        CALL SHOWHLP('detailed description')
        LSIGCR=.FALSE.
        LFACTCR=.FALSE.
        LMINVCR=.FALSE.
        LOTHERCR=.FALSE.
        NFIT=2
        NDP=2
        NPP=25
        PERCENT1=10.
        PERCENT2=60.
        RQSKY=60.
        NTOTCR=0
        NCRSAVED=0
        FACTOR=100.
        NAUXFRAME=0
        ANGLE3D=0
C
        WRITE(*,100)'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
C
        WRITE(*,100)'Input file name'
        NAMFIL=INFILEX(20,'@',NS,NC,STWV,DISP,1,.FALSE.)
        WRITE(*,100)'Reading file...'
        DO I=1,NS
          READ(20,ERR=910)(A(J,I),J=1,NC)
        END DO
        CLOSE(20)
        L1=TRUELEN(NAMFIL)
        L2=TRUELEN(OBJECT)
        WRITE(*,101)'  ...OK!'
        IF(CERR.EQ.'y')THEN
          WRITE(*,100)'Error file name '
          CALL GUESSEF(NAMFIL,ERRFILE)
          ERRFILE=INFILEX(30,ERRFILE,NS,NC,STWV,DISP,21,.TRUE.)    !match
          WRITE(*,100)'Reading file...'
          DO I=1,NS
            READ(30,ERR=910) (ERR(J,I),J=1,NC)
          END DO
          CLOSE(30)
          WRITE(*,101)'  ...OK!'
        END IF
        IF(L2.GT.0) NAMFIL(L1+1:L1+L2+3)=' ['//OBJECT(1:L2)//']'
C duplicate variable names (to be compatible with other Subroutines)
        NSCAN=NS
        NCHAN=NC
        NC1=1
        NC2=NC
        NS1=1
        NS2=NS
C
        WRITE(*,*)
        WRITE(*,101)'* DEFAULT OPTIONS:'
        WRITE(*,101)'SIGCR = 5.0'
        WRITE(*,101)'SIGTHRESHOLD = 6.0'
        WRITE(*,101)'FACTCR (not employed)'
        WRITE(*,101)'MINVCR (not employed)'
        WRITE(*,101)'SIGCRAUX = 10.0'
        WRITE(*,101)'Do remove mean x/y direction before '//
     +   'looking for C.R.'
        WRITE(*,101)'Edge size of the Searching square: 15 pixels'
        WRITE(*,101)'Plot type 1'
        WRITE(*,101)'FG and BG do not fixed'
        CALL SHOWHLP('default options')
        WRITE(*,101)
        WRITE(*,100)'Accept default options (y/n) '
        CDEF(1:1)=READC('y','yn')
        LSIGCR=.TRUE.
        LFACTCR=.FALSE.
        LMINVCR=.FALSE.
        LOTHERCR=.FALSE.
        LREMOVEX=.TRUE.
        LREMOVEY=.TRUE.
        SIGCR=5.0
        SIGTHRESHOLD=6.0
        SIGCRAUX=10.0
        NWIN=15
        NTYPEPLOT=1
        CFIX_FGBG='n'
        IF(CDEF.EQ.'n')THEN
          CALL CROPTIONS
        END IF
C
        CALL SHOWHLP('history file')
        WRITE(*,100)'Create history file... (y/n) '
        CDEF(1:1)=READC('n','yn')
        IF(CDEF.EQ.'y')THEN
          LHISTFILE=.TRUE.
          WRITE(*,100)'History file name'
          HISTFILE=OUTFILEX(80,'@',0,0,0.,0.,3,.FALSE.)
        ELSE
          LHISTFILE=.FALSE.
        END IF
C
        CALL RPGBEGIN(NTERM,IDN,LCOLOR)
        CALL PALETTE(3)
        CALL BUTTSYB(3)
        CALL BUTTSIT(.TRUE.)
C
14      CALL PGBBUF
        CALL PGSAVE
        DO I=1,NBOTONES
          CALL BUTTON(I,LABEL(I),BMODE(I))
        END DO
        CALL PGEBUF
        CALL PGUNSA
C
15      WRITE(*,101)
        WRITE(*,101)'---------------------  MAIN MENU   --------'//
     +   '------------------'
        WRITE(*,101)'start.....- begin automatic detection of C.R. '//
     +   '(clean by HAND)'
        WRITE(*,101)'region....- examination of some pixel region'
        WRITE(*,101)'window....- change display window edge size'
        WRITE(*,101)'automatic.- clean automatically'
        WRITE(*,101)'look......- have a look to the image'
        WRITE(*,101)'options...- change searching options'
        WRITE(*,101)'save......- save current image'
        WRITE(*,101)'histogram1- create DATA histogram'
        WRITE(*,101)'histogram2- create SIGMA histogram'
        WRITE(*,101)'top100....- search the top 100 in SIGMA'
        WRITE(*,101)'plotsp3d  - emulate plotsp3d program'
        WRITE(*,101)'QUIT......- end of program'
        WRITE(*,101)'-------------------------------------------'//
     +   '------------------'
        WRITE(*,101) ' '
        WRITE(*,101) 'NOTE: remember that # reverses data!!!'
        WRITE(*,101) ' '
        CALL SHOWHLP('menu principal')
16      CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        IF(CH.EQ.'s')THEN
          IMENU=1
        ELSEIF(CH.EQ.'r')THEN
          IMENU=2
        ELSEIF(CH.EQ.'w')THEN
          IMENU=3
        ELSEIF(CH.EQ.'a')THEN
          IMENU=4
        ELSEIF(CH.EQ.'l')THEN
          IMENU=5
        ELSEIF(CH.EQ.'o')THEN
          IMENU=6
        ELSEIF(CH.EQ.'v')THEN
          IMENU=0
        ELSEIF(CH.EQ.'1')THEN
          IMENU=8
        ELSEIF(CH.EQ.'2')THEN
          IMENU=9
        ELSEIF(CH.EQ.'t')THEN
          IMENU=10
        ELSEIF(CH.EQ.'d')THEN
          IMENU=11
        ELSEIF(CH.EQ.'q')THEN
          IMENU=-1
        ELSEIF(CH.EQ.'#')THEN
          WRITE(*,100) 'Reversing data...'
          DO I=1,NS
            DO J=1,NC
              A(J,I)=-A(J,I)
            END DO
          END DO
          WRITE(*,*) '...OK!'
          GOTO 16
        ELSE
          IF(NB.EQ.0)GOTO 16
          IMENU=NB
          IF(NB.EQ.12) IMENU=-1
          IF(NB.EQ.7) IMENU=0
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.-1)THEN
          CALL BUTTON(12,LABEL(12),5) 
          IF(NTOTCR.GT.NCRSAVED)THEN
            WRITE(*,101)'WARNING: last changes have not been saved.'
            WRITE(*,100)'Do you really want to QUIT (y/n) '
            CSURE(1:1)=READC('n','yn')
            IF(CSURE.NE.'y')THEN
              CALL BUTTON(12,LABEL(12),0) 
              GOTO 15
            END IF
          END IF
          CALL PGEND
          IF(LHISTFILE)THEN
            WRITE(80,'(A,I5)')'No. of C.R. removed: ',NTOTCR
            CLOSE(80)
          END IF
          WRITE(*,'(A,I5)')'No. of C.R. removed: ',NTOTCR
          STOP
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.0)THEN
          CALL BUTTON(7,LABEL(7),5)
          WRITE(*,100)'* Previous name: '
          WRITE(*,101) NAMFIL(1:TRUELEN(NAMFIL))
          WRITE(*,100)'Output file name'
          NAMFIL=OUTFILEX(30,'@',NS,NC,STWV,DISP,1,.FALSE.)
          L1=TRUELEN(NAMFIL)
          L2=TRUELEN(OBJECT)
          WRITE(*,100)'Saving file...'
          DO I=1,NS
            WRITE(30)(A(J,I),J=1,NC)
          END DO
          CLOSE(30)
          WRITE(*,101)'...OK! (file saved and closed)'
          IF(CERR.EQ.'y')THEN
            WRITE(*,100)'Output error file name '
            CALL GUESSEF(NAMFIL,ERRFILE)
            ERRFILE=OUTFILEX(32,ERRFILE,NS,NC,STWV,DISP,1,.TRUE.)
            WRITE(*,100)'Saving file...'
            DO I=1,NS
              WRITE(32)(ERR(J,I),J=1,NC)
            END DO
            CLOSE(32)
            WRITE(*,101)'...OK! (file saved and closed)'
          END IF
          IF(L2.GT.0) NAMFIL(L1+1:L1+L2+3)=' ['//OBJECT(1:L2)//']'
          WRITE(*,*)
          NCRSAVED=NTOTCR
          CALL BUTTON(7,LABEL(7),0)
          GOTO 15
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.3)THEN
          CALL BUTTON(3,LABEL(3),5)
          WRITE(*,101)'New size?'
          WRITE(*,100)'(Note: this number must be odd) '
          WRITE(CDUMMY,*)NPP
          NPP=READI(CDUMMY)
          IF(MOD(NPP,2).EQ.0)THEN
            NPP=NPP+1
            WRITE(*,100)'WARNING: Effective Edge Size changed to '
            WRITE(*,*)NPP
          END IF
          IF(NPP.GT.NPPMAX)THEN
            WRITE(*,100)'NPP value set to NPPMAX: '
            WRITE(*,*)NPPMAX
            NPP=NPPMAX
          END IF
          CALL BUTTON(3,LABEL(3),0)
          GOTO 15
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.6)THEN
          CALL BUTTON(6,LABEL(6),5)
          CALL CROPTIONS
          CALL BUTTON(6,LABEL(6),0)
          GOTO 15
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.5)THEN
          CLEANAUTO=.FALSE.
          CALL BUTTON(5,LABEL(5),5)
          CALL PGERAS
ccc       CALL RPGERAS
ccc       CALL RPGERASB
          CALL BUTTSPR(0.10,0.95,0.10,0.70)
C
          CALL LOOK(NAMFIL,NC1,NC2,NS1,NS2)
C
          CALL RPGERASB
          GOTO 14
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.8)THEN
          CALL BUTTON(8,LABEL(8),5)
          CALL HISTOGRAMA(1,NAMFIL)
          CALL BUTTON(8,LABEL(8),0)
          GOTO 15
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.9)THEN
          CALL BUTTON(9,LABEL(9),5)
          CALL HISTOGRAMA(2,NAMFIL)
          CALL BUTTON(9,LABEL(9),0)
          GOTO 15
        END IF
C------------------------------------------------------------------------------
C Nota: la llamada a HISTOGRAMA tambien inicializa las variables SCANYES,
C CHANYES, y TOP100CLEAN, de modo que podemos alterar despues estos valores
C en SUBPLOT.
        IF(IMENU.EQ.10)THEN
          CALL BUTTON(10,LABEL(10),5)
          CALL HISTOGRAMA(3,NAMFIL)
          CALL BUTTON(10,LABEL(10),0)
C
          CALL PGBBUF
          CALL PGSAVE
          CALL PGERAS
          DO K=1,NBOTONES
            CDUMMY=LABEL2(K)
            IF(K.EQ.6)THEN
              WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL BUTTON(K,CDUMMY(1:L),0)
            ELSE
              IF(BMODE2(K).GE.0) CALL BUTTON(K,CDUMMY,0)
              IF(BMODE2(K).GT.0) CALL BUTTON(K,CDUMMY,BMODE2(K))
            END IF
          END DO
          CALL PGEBUF
          CALL PGUNSA
C
          WRITE(*,*)
          LNEXT=.TRUE.
          DO NTOP=NTOPMAX,1,-1
            IF(LNEXT)THEN
              I=ITOP(NTOP)
              J=JTOP(NTOP)
            ELSE
              I=0
              J=0
            END IF
            IF((I.NE.0).AND.(J.NE.0))THEN
              WRITE(*,'(A,I3.3,A)')'>>> Top#',NTOP,':'
              WRITE(CDUMMY,'(I10,A1,I10)')J,',',I
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              WRITE(*,100)'channel, scan: '
              WRITE(*,101)CDUMMY(1:L)
              WRITE(*,100)'ATOP (in sigma units): '
              WRITE(*,*) ATOP(NTOP)
              IF(SCANYES(I).AND.CHANYES(J))THEN
                IF(TOP100CLEAN(J,I))THEN
                  WRITE(*,100)'WARNING: pixel skipped'
                  WRITE(*,101)' (already cleaned).'
                ELSE
                  CALL STATISTICS(J,I)
                  CALL RPGERASW(0.,1.,0.,0.80)
                  CALL SUBPLOT(J,I,LNEXT,.TRUE.)
                END IF
              ELSE
                WRITE(*,101)'WARNING: pixel skipped.'
              END IF
            END IF
          END DO
C
          CALL RPGERASB
          GOTO 14
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.11)THEN
          CALL BUTTON(10,LABEL(10),5)
          CALL RPGERASB
          CALL BUTTSYB(2)
          CALL PLOTSP3D(NAMFIL)
          CALL RPGERASB
          CALL BUTTSYB(3)
          GOTO 14
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.2)THEN
          CALL BUTTON(2,LABEL(2),5)
17        WRITE(*,100)'Channel, scan to go'
          CALL READ2I('@',J,I)
          IF((I.LT.1).OR.(I.GT.NS).OR.(J.LT.1).OR.(J.GT.NC))THEN
            WRITE(*,101)'ERROR: pixel outside Image. Try again.'
            GOTO 17
          END IF
          NWIN=NPP
          CALL STATISTICS(J,I)
          CALL PGBBUF
          CALL PGSAVE
          CALL PGERAS
ccc       CALL RPGERAS
ccc       CALL RPGERASB
ccc       CALL PGIDEN_RED
          DO K=1,NBOTONES
            CDUMMY=LABEL2(K)
            IF(K.EQ.6)THEN
              WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL BUTTON(K,CDUMMY(1:L),0)
            ELSE
              IF(BMODE2(K).GE.0) CALL BUTTON(K,CDUMMY,0)
              IF(BMODE2(K).GT.0) CALL BUTTON(K,CDUMMY,BMODE2(K))
            END IF
          END DO
          CALL PGEBUF
          CALL PGUNSA
          CALL SUBPLOT(J,I,LNEXT,.TRUE.)
          CALL RPGERASB
          GOTO 14
        END IF
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C                          IMENU=1 or IMENU=4
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        CALL BUTTON(IMENU,LABEL(IMENU),5)
C SCANYES(I) determina si en el scan numero I se procedera o no
C a la busqueda automatica de rayos cosmicos para IMENU=1 o 4.
C CHANYES(J,I) determina si en el canal J se procedera a la busqueda de
C rayos cosmicos.
C..............................................................................
        NC1P=1
        NC2P=NC
        NS1P=1
        NS2P=NS
        CALL BUTTSPR(0.10,0.95,0.10,0.70)
20      CALL RPGERASW(0.,1.,0.,0.80)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGENV(REAL(NC1P)-.6,REAL(NC2P)+.6,REAL(NS1P)-.6,
     +      REAL(NS2P)+.6,0,0)
          ELSE
            CALL PGENV(REAL(NC1P)-.6,REAL(NC2P)+.6,REAL(NS1P)-.6,
     +      REAL(NS2P)+.6,0,-2)
ccc            CALL PGWINDOW(REAL(NC1P)-.6,REAL(NC2P)+.6,REAL(NS1P)-.6,
ccc     +      REAL(NS2P)+.6)
            CALL PGBOX('BCNITS',0.0,0,'BCNITS',0.0,0)
          END IF
          CALL PGIDEN_RED
        END DO
        CALL STATISTICS((NC1P+NC2P)/2,(NS1P+NS2P)/2)
        IF(LSIGCR)THEN
          BG=MEAN-SIGCR*SIGMA
          FG=MEAN+SIGCR*SIGMA
        ELSEIF(LFACTCR)THEN
          BG=(MEAN-1)*FACTCR
          FG=(MEAN+1)*FACTCR
        ELSEIF(LMINVCR)THEN
          BG=MEAN+MINVCR
          FG=MEAN-MINVCR
        ELSEIF(LOTHERCR)THEN
          BG=MEAN-SIGCRAUX*SIGMA
          FG=MEAN+SIGCRAUX*SIGMA
        ELSE
          STOP 'FATAL ERROR: you should never reach this point!'
        END IF
        MINVAL=A(NC1P,NS1P)
        MAXVAL=MINVAL
        DO I=NS1P,NS2P
          DO J=NC1P,NC2P
            IF(A(J,I).LT.MINVAL) MINVAL=A(J,I)
            IF(A(J,I).GT.MAXVAL) MAXVAL=A(J,I)
          END DO
        END DO
!       IF(BG.LT.MINVAL) BG=MINVAL
!       IF(FG.GT.MAXVAL) FG=MAXVAL
        CALL ZSCALE(A,NC1,NC2,NS1,NS2,BG,FG)
        IF(BG.EQ.FG)THEN
          IF(BG.EQ.0.)THEN
            BG=-1.
            FG=+1.
          ELSE
            BG=0.9*BG
            FG=1.1*FG
          END IF
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGIMAG(A,NCMAX,NSMAX,NC1P,NC2P,NS1P,NS2P,FG,BG,TR)
          CALL PGMTEXT('L',2.5,.5,.5,'scan')
          CALL PGMTEXT('B',3.,.5,.5,'channel')
          CALL PGMTEXT('T',1.,0.,0.,'file: '//NAMFIL)
        END DO
C
C..............................................................................
25      WRITE(*,101)'(1) Change plot limits'
        WRITE(*,101)'(2) Change BG/FG'
        WRITE(*,101)'(3) Select plotted image and remove '//
     +   'scans/channels with mouse'
        WRITE(*,101)'(4) Select regions to be cleaned with keyboard'
        WRITE(*,100)'Option '
        CMOUSE(1:1)=READC('4','1234')
C..............................................................................
        IF(CMOUSE.EQ.'1')THEN
          WRITE(*,101)'* Select image region:'
          WRITE(CDUMMY,110)'1,',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'First and last scan '
          CALL READ2I(CDUMMY(1:L),NS1P,NS2P)
          WRITE(CDUMMY,110)'1,',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,100)'First and last channel '
          CALL READ2I(CDUMMY(1:L),NC1P,NC2P)
          GOTO 20
C..............................................................................
        ELSEIF(CMOUSE.EQ.'2')THEN
          WRITE(CDUMMY,*)BG
          WRITE(*,100)'Background '
          BG=READF(CDUMMY)
          WRITE(CDUMMY,*)FG
          WRITE(*,100)'Foreground '
          FG=READF(CDUMMY)
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGIMAG(A,NCMAX,NSMAX,NC1P,NC2P,NS1P,NS2P,FG,BG,TR)
          END DO
          GOTO 25
C..............................................................................
        ELSEIF(CMOUSE.EQ.'3')THEN
          DO I=1,NS
            SCANYES(I)=.FALSE.
          END DO
          DO I=NS1P,NS2P
            SCANYES(I)=.TRUE.               !por definicion todos los dibujados
          END DO
          DO J=1,NC
            CHANYES(J)=.FALSE.
          END DO
          DO J=NC1P,NC2P
            CHANYES(J)=.TRUE.               !por definicion todos los dibujados
          END DO
          CH=' '
          DO WHILE(CH.NE.'X')
            WRITE(*,100)'Remove scans with mouse...'
            IF(LCOLOR(1)) CALL PGSCI(2)
            CALL RPGBAND(5,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            IF(CH.NE.'X')THEN
              NS1YES=NINT(YC)
              IF(NS1YES.LT.NS1P) NS1YES=NS1P
              IF(NS1YES.GT.NS2P) NS1YES=NS2P
              WRITE(*,100)'  press again...'
              IF(LCOLOR(1)) CALL PGSCI(2)
              CALL RPGBAND(3,0,0.,REAL(NS1YES),XC,YC,CH)
              IF(LCOLOR(1)) CALL PGSCI(1)
              WRITE(*,*)
              IF(CH.NE.'X')THEN
                NS2YES=NINT(YC)
                IF(NS2YES.LT.NS1P) NS2YES=NS1P
                IF(NS2YES.GT.NS2P) NS2YES=NS2P
                IF(NS2YES.LT.NS1YES)THEN
                  NS0YES=NS2YES
                  NS2YES=NS1YES
                  NS1YES=NS0YES
                END IF
                DO I=NS1YES,NS2YES
                  SCANYES(I)=.FALSE.
                END DO
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM))CALL PGSCI(5)
                  CALL PGMOVE(REAL(NC1P),REAL(NS1YES))
                  CALL PGDRAW(REAL(NC2P),REAL(NS1YES))
                  CALL PGMOVE(REAL(NC1P),REAL(NS2YES))
                  CALL PGDRAW(REAL(NC2P),REAL(NS2YES))
                  CALL PGRECT(REAL(NC1P),REAL(NC1P)+
     +             REAL(NC2P-NC1P+1)/100.,REAL(NS1YES),REAL(NS2YES))
                  CALL PGRECT(REAL(NC2P),REAL(NC2P)-
     +             REAL(NC2P-NC1P+1)/100.,REAL(NS1YES),REAL(NS2YES))
                  IF(LCOLOR(ITERM))CALL PGSCI(1)
                END DO
              ELSE
                WRITE(*,*)
              END IF
            END IF
          END DO
          WRITE(*,*)
          CH=' '
          DO WHILE(CH.NE.'X')
            WRITE(*,100)'Remove channels with mouse...'
            IF(LCOLOR(1)) CALL PGSCI(2)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            IF(CH.NE.'X')THEN
              NC1YES=NINT(XC)
              IF(NC1YES.LT.NC1P) NC1YES=NC1P
              IF(NC1YES.GT.NC2P) NC1YES=NC2P
              WRITE(*,100)'  press again...'
              IF(LCOLOR(1)) CALL PGSCI(1)
              CALL RPGBAND(4,0,REAL(NC1YES),0.,XC,YC,CH)
              IF(LCOLOR(1)) CALL PGSCI(2)
              WRITE(*,*)
              IF(CH.NE.'X')THEN
                NC2YES=NINT(XC)
                IF(NC2YES.LT.NC1P) NC2YES=NC1P
                IF(NC2YES.GT.NC2P) NC2YES=NC2P
                IF(NC2YES.LT.NC1YES)THEN
                  NC0YES=NC2YES
                  NC2YES=NC1YES
                  NC1YES=NC0YES
                END IF
                DO J=NC1YES,NC2YES
                  CHANYES(J)=.FALSE.
                END DO
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  IF(LCOLOR(ITERM))CALL PGSCI(4)
                  CALL PGMOVE(REAL(NC1YES),REAL(NS1P))
                  CALL PGDRAW(REAL(NC1YES),REAL(NS2P))
                  CALL PGMOVE(REAL(NC2YES),REAL(NS1P))
                  CALL PGDRAW(REAL(NC2YES),REAL(NS2P))
                  CALL PGRECT(REAL(NC1YES),REAL(NC2YES),
     +             REAL(NS1P),REAL(NS1P)+REAL(NS2P-NS1P+1)/100.)
                  CALL PGRECT(REAL(NC1YES),REAL(NC2YES),
     +             REAL(NS2P),REAL(NS2P)-REAL(NS2P-NS1P+1)/100.)
                  IF(LCOLOR(ITERM))CALL PGSCI(1)
                END DO
              ELSE
                WRITE(*,*)
              END IF
            END IF
          END DO
          WRITE(*,*)
C..............................................................................
        ELSEIF(CMOUSE.EQ.'4')THEN
          DO I=1,NS
            SCANYES(I)=.FALSE.                          !por definicion ninguno
          END DO
          WRITE(*,110)'* Define scans to be cleaned: '//
     +     'valid range from 1 to ',NSCAN
          LEXIT=.FALSE.
          DO WHILE(.NOT.LEXIT)
            WRITE(*,100)'Scan region (0,0=EXIT) '
            CALL READ2I('0,0',NS1YES,NS2YES)
            IF((NS1YES.EQ.0).AND.(NS2YES.EQ.0))THEN
              DO I=1,NS
                IF(SCANYES(I)) LEXIT=.TRUE.
              END DO
              IF(.NOT.LEXIT)THEN
                WRITE(*,101)'ERROR: no scan region has been defined.'
              END IF
            ELSE
              IF(.NOT.INSIDE2(NS1YES,NS2YES,1,NS))THEN
                WRITE(*,101)'ERROR: Invalid scan region. Try again.'
              ELSE
                DO I=NS1YES,NS2YES
                  SCANYES(I)=.TRUE.
                END DO
              END IF
            END IF
          END DO
          DO J=1,NC
            CHANYES(J)=.FALSE.                          !por definicion ninguno
          END DO
          WRITE(*,110)'* Define channels to be cleaned: '//
     +     'valid range from 1 to ',NCHAN
          LEXIT=.FALSE.
          DO WHILE(.NOT.LEXIT)
            WRITE(*,100)'Channel region (0,0=EXIT) '
            CALL READ2I('0,0',NC1YES,NC2YES)
            IF((NC1YES.EQ.0).AND.(NC2YES.EQ.0))THEN
              DO J=1,NC
                IF(CHANYES(J)) LEXIT=.TRUE.
              END DO
              IF(.NOT.LEXIT)THEN
                WRITE(*,101)'ERROR: no channel region has been defined.'
              END IF
            ELSE
              IF(.NOT.INSIDE2(NC1YES,NC2YES,1,NC))THEN
                WRITE(*,101)'ERROR: Invalid channel region. Try again.'
              ELSE
                DO J=NC1YES,NC2YES
                  CHANYES(J)=.TRUE.
                END DO
              END IF
            END IF
          END DO
C..............................................................................
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100)'Skip channels by keyboard (y/n) '
        CSKIPCH(1:1)=READC('n','yn')
        IF(CSKIPCH.EQ.'y')THEN
          WRITE(*,110)'* Valid range from 1 to ',NCHAN
          LEXIT=.FALSE.
          DO WHILE(.NOT.LEXIT)
            WRITE(*,100)'Channel region (0,0=EXIT) '
            CALL READ2I('0,0',NC1YES,NC2YES)
            IF((NC1YES.EQ.0).AND.(NC2YES.EQ.0))THEN
              LEXIT=.TRUE.
            ELSE
              IF(.NOT.INSIDE2(NC1YES,NC2YES,1,NC))THEN
                WRITE(*,101)'ERROR: Invalid channel region. Try again.'
              ELSE
                DO J=NC1YES,NC2YES
                  CHANYES(J)=.FALSE.
                END DO
              END IF
            END IF
          END DO
        END IF
C------------------------------------------------------------------------------
        IF(IMENU.EQ.4)THEN
          WRITE(*,101)'(1) X interpolation'
          WRITE(*,101)'(2) Y interpolation'
          WRITE(*,101)'(3) Polynomial surface'
          IF(NAUXFRAME.GT.0)THEN
            WRITE(*,101)'(4) Replace by auxiliary frame(s)'
          END IF
          WRITE(*,101)'(0) NONE (RETURN to main menu)'
          WRITE(*,100)'Option '
          CAUTOMODE=READILIM('2',0,4)
          IF(CAUTOMODE.EQ.0)THEN
            CALL BUTTON(4,LABEL(4),0)
            GOTO 15
          END IF
          CLEANAUTO=.TRUE.
          CCROS='n'
          WRITE(*,100)'Plot individual c.r. (y/n) '
          CPINTAAUTO(1:1)=READC('n','yn')
          IF(CPINTAAUTO.EQ.'n')THEN
            WRITE(*,100)'Plot crosses over c.r. (y/n) '
            CCROS(1:1)=READC('y','yn')
          END IF
        ELSE
          CLEANAUTO=.FALSE.
          CCROS='n'
        END IF
C------------------------------------------------------------------------------
        IF(CCROS.EQ.'n')THEN
          CALL BUTTON(IMENU,LABEL(IMENU),0)
          CALL PGERAS
ccc       CALL RPGERAS
ccc       CALL RPGERASB
          DO I=1,NBOTONES
            CDUMMY=LABEL2(I)
            IF(I.EQ.6)THEN
              WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL BUTTON(I,CDUMMY(1:L),0)
            ELSE
              IF(BMODE2(I).GE.0) CALL BUTTON(I,CDUMMY,0)
              IF(BMODE2(I).GT.0) CALL BUTTON(I,CDUMMY,BMODE2(I))
            END IF
          END DO
        END IF
C
33      DO I=1,NS
          IF(SCANYES(I))THEN
            !WRITE(*,'(A,I5.5,$)')'\\b\\b\\b\\b\\b',I
            DO J=1,NC
              IF(CHANYES(J))THEN
                CRFOUND=.FALSE.
                CALL STATISTICS(J,I)
                IF(LSIGCR)THEN
                  DEV=(A(J,I)-MEAN)/SIGMA
                  IF(DEV.GT.SIGCR) CRFOUND=.TRUE.
                END IF
                IF(LFACTCR)THEN
                  DEV=(A(J,I)-MEAN)/MEAN
                  IF(DEV.GT.FACTCR) CRFOUND=.TRUE.
                END IF
                IF(LMINVCR)THEN
                  DEV=A(J,I)-MEAN
                  IF(DEV.GT.MINVCR) CRFOUND=.TRUE.
                END IF
                IF(LOTHERCR)THEN
                  CALL CHECK_AUXFRAME(J,I,DEV)
                  IF(DEV.GT.SIGCRAUX) CRFOUND=.TRUE.
                END IF
                IF(CRFOUND)THEN
                  WRITE(*,*)
                  WRITE(*,105)'C.R. FOUND AT ',J,I
                  WRITE(*,100)'Pixel value: '
                  WRITE(*,*)A(J,I)
                  WRITE(*,100)'MEAN,SIGMA : '
                  WRITE(*,*)MEAN,SIGMA
                  WRITE(*,100)'SIGMA times over MEAN: '
                  WRITE(*,*)(A(J,I)-MEAN)/SIGMA
                  IF(CLEANAUTO)THEN
                    IF(CPINTAAUTO.EQ.'n')THEN
                      IF(CCROS.EQ.'y')THEN
                        DO ITERM=NTERM,1,-1
                          CALL PGSLCT(IDN(ITERM))
                          IF(LCOLOR(ITERM))CALL PGSCI(3)
                          CALL PGPOINT(1,REAL(J),REAL(I),23)
                          IF(LCOLOR(ITERM))CALL PGSCI(1)
                        END DO
                      END IF
                    END IF
                  END IF
                  IF(CCROS.EQ.'n')THEN
                    CALL RPGERASW(0.,1.,0.,0.80)
ccc                 CALL RPGERAS
ccc                 CALL PGIDEN_RED
                    CALL SUBPLOT(J,I,LNEXT,.TRUE.)
                  ELSE
                    CALL SUBPLOT(J,I,LNEXT,.FALSE.)
                  END IF
                  IF(.NOT.SCANYES(I)) GOTO 34
                  IF((.NOT.LNEXT).AND.(IMENU.EQ.1))THEN
                    CALL RPGERASB
                    GOTO 14
                  END IF
                END IF
              END IF                                          !IF en CHANYES(J)
            END DO                                                 !bucle en NC
34          CONTINUE
          END IF                                              !IF en SCANYES(I)
        END DO                                                     !bucle en NS
C
        WRITE(*,*)
        IF((IMENU.EQ.4).AND.(CMOUSE.EQ.'3'))THEN
          WRITE(*,100)'Repeat cleaning in the same regions (y/n) '
          CCLEAN2(1:1)=READC('n','yn')
          IF(CCLEAN2.EQ.'y') GOTO 33
        END IF
        CLEANAUTO=.FALSE.
C
        IF(CCROS.EQ.'n')THEN
          CALL RPGERASB
          GOTO 14
        ELSE
          IF(IMENU.EQ.4) CALL BUTTON(IMENU,LABEL(IMENU),0)
          GOTO 15
        END IF
C
105     FORMAT(A,I5,2X,I5)
C
910     WRITE(*,101)'ERROR reading file '//NAMFIL
        CLOSE(13)
        STOP
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE STATISTICS(J,I)
C Esta subrutina calcula la media y la desviacion standard en el cuadrado
C de busqueda teniendo en cuenta el efecto del borde. Asimismo, en este
C calculo no se tiene en cuenta el SCAN ni el CANAL donde se encuentra
C el pixel (J,I) que estamos examinando. El retorno se realiza a traves de
C las variables MEAN y SIGMA.
C LREMOVEX=.TRUE.: elimina forma en direccion x antes de calcular MEAN,SIGMA
C LREMOVEY=.TRUE.: elimina forma en direccion y antes de calcular MEAN,SIGMA
        IMPLICIT NONE
        INTEGER I,J
C
        INCLUDE 'redlib.inc'
C
        INTEGER NS,NC,NWIN
        INTEGER NS1,NS2,NC1,NC2
        INTEGER NWIN2
        INTEGER II,JJ,NPIX
        REAL MEAN,SIGMA
        REAL A(NCMAX,NSMAX)
        REAL MEANX(NCMAX),MEANY(NSMAX)
        LOGICAL LREMOVEX,LREMOVEY
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLK1/MEAN,SIGMA
        COMMON/BLK2/NS,NC,NWIN
        COMMON/BLK18/LREMOVEX,LREMOVEY
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        NWIN2=(NWIN-1)/2
        NS1=I-NWIN2
        NS2=I+NWIN2
        NC1=J-NWIN2
        NC2=J+NWIN2
        IF(NS1.LT.1)THEN
          NS1=1
          NS2=NWIN
        END IF
        IF (NS2.GT.NS)THEN
          NS1=NS-NWIN+1
          NS2=NS
        END IF
        IF (NC1.LT.1)THEN
          NC1=1
          NC2=NWIN
        END IF
        IF (NC2.GT.NC)THEN
          NC1=NC-NWIN+1
          NC2=NC
        END IF
C
        NPIX=(NWIN-1)*(NWIN-1)
C
C forma promedio en x
        DO JJ=NC1,NC2
          MEANX(JJ)=0.
        END DO
        IF(LREMOVEX)THEN
          DO JJ=NC1,NC2
            DO II=NS1,NS2
              IF(II.NE.I)THEN
                MEANX(JJ)=MEANX(JJ)+A(JJ,II)
              END IF
            END DO
            MEANX(JJ)=MEANX(JJ)/REAL(NWIN-1)
          END DO
        END IF
C forma promedio en y
        DO II=NS1,NS2
          MEANY(II)=0.
        END DO
        IF(LREMOVEY)THEN
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(JJ.NE.J)THEN
                MEANY(II)=MEANY(II)+A(JJ,II)-MEANX(JJ)
              END IF
            END DO
            MEANY(II)=MEANY(II)/REAL(NWIN-1)
          END DO
        END IF
C
        MEAN=0.
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            IF((II.NE.I).AND.(JJ.NE.J))THEN
              MEAN=MEAN+(A(JJ,II)-MEANX(JJ)-MEANY(II))
            END IF
          END DO
        END DO
        MEAN=MEAN/REAL(NPIX)
C
        SIGMA=0.
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            IF((II.NE.I).AND.(JJ.NE.J))THEN
              SIGMA=SIGMA+(MEAN-A(JJ,II)+MEANX(JJ)+MEANY(II))*
     +         (MEAN-A(JJ,II)+MEANX(JJ)+MEANY(II))
            END IF
          END DO
        END DO
        SIGMA=SQRT(SIGMA/REAL(NPIX-1))
C
        MEAN=MEAN+MEANX(J)+MEANY(I)                  !para que la cosa funcione
C
        END
C
C******************************************************************************
C
        SUBROUTINE HISTOGRAMA(NTYPE,NAMFIL)
C Realiza un histograma de la imagen. Segun el valor de NTYPE:
C 1 - Histograma en senhal
C 2 - Histograma en desviaciones respecto a la media en unidades de SIGMA local
C 3 - Top 100: busca los 100 pixels con mayor valor en la imagen de SIGMA local
        IMPLICIT NONE
        INTEGER NTYPE
        CHARACTER*(*) NAMFIL
C
        INTEGER NBINMAX
        PARAMETER (NBINMAX=1000)                     !numero maximo de binnings
        INTEGER NTOPMAX
        PARAMETER (NTOPMAX=1000)                     !numero maximo para TOP100
C
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J,K,L
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NY(NBINMAX)
        INTEGER NBIN
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NTOP,ITOP(NTOPMAX),JTOP(NTOPMAX)
        REAL MEAN,SIGMA
        REAL A(NCMAX,NSMAX)
        REAL B(NCMAX,NSMAX)
        REAL BMIN,BMAX,DB,BMINF,BMAXF
        REAL X(NBINMAX),Y(NBINMAX)
        REAL YMIN,YMAX,DY
        REAL ATOP(NTOPMAX)
        CHARACTER*1 CREPLOT
        CHARACTER*75 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL SCANYES(NSMAX),CHANYES(NCMAX)
        LOGICAL TOP100CLEAN(NCMAX,NSMAX)
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLK1/MEAN,SIGMA
        COMMON/BLK19/CHANYES,SCANYES
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
        COMMON/BLK1TOP100/ITOP,JTOP
        COMMON/BLK2TOP100/ATOP
        COMMON/BLK3TOP100/TOP100CLEAN
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C seleccionamos region
        IF(NTYPE.EQ.3)THEN
          WRITE(*,101)'Enter region to be used to calculate top 100:'
        ELSE
          WRITE(*,101)'Enter region to be used to calculate histogram:'
        END IF
10      WRITE(CDUMMY,'(A2,I10)') '1,',NSCAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'1st and last scan '
        CALL READ2I(CDUMMY(1:L),NS1,NS2)
        IF((NS1.LT.1).OR.(NS2.GT.NSCAN).OR.(NS1.GT.NS2))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 10
        END IF
12      WRITE(CDUMMY,'(A2,I10)') '1,',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'1st and last channel '
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        IF((NC1.LT.1).OR.(NC2.GT.NCHAN).OR.(NC1.GT.NC2))THEN
          WRITE(*,101)'ERROR: numbers out of range. Try again.'
          GOTO 12
        END IF
C
        IF(NTYPE.EQ.3)THEN
          DO I=1,NSCAN
            SCANYES(I)=.FALSE.
          END DO
          DO I=NS1,NS2
            SCANYES(I)=.TRUE.
          END DO
          DO J=1,NCHAN
            CHANYES(J)=.FALSE.
          END DO
          DO J=NC1,NC2
            CHANYES(J)=.TRUE.
          END DO
          DO I=1,NSCAN
            DO J=1,NCHAN
              TOP100CLEAN(J,I)=.FALSE.
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
C seleccionamos datos a representar en el histograma
        WRITE(*,100)'Wait...'
        IF(NTYPE.EQ.1)THEN
          DO I=NS1,NS2
            DO J=NC1,NC2
              B(J,I)=A(J,I)
            END DO
          END DO
        ELSEIF((NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
          WRITE(*,*)
          DO I=NS1,NS2
            !WRITE(*,'(A,I4.4,$)')'\\b\\b\\b\\b',I
            DO J=NC1,NC2
              CALL STATISTICS(J,I)
              B(J,I)=(A(J,I)-MEAN)/SIGMA
            END DO
          END DO
          WRITE(*,*)
        END IF
        WRITE(*,101)'  ..OK!'
C------------------------------------------------------------------------------
C calculamos TOP100
        IF(NTYPE.EQ.3)THEN
          WRITE(*,101)'Sorting...'
          DO NTOP=1,NTOPMAX
            ITOP(NTOP)=0
            JTOP(NTOP)=0
          END DO
          NTOP=0
          DO I=NS1,NS2
            !WRITE(*,'(A,I4.4,$)')'\\b\\b\\b\\b',I
            DO J=NC1,NC2
              IF(NTOP.LT.NTOPMAX)THEN
                NTOP=NTOP+1
                ITOP(NTOP)=I
                JTOP(NTOP)=J
                ATOP(NTOP)=B(J,I)
                IF(NTOP.EQ.NTOPMAX)THEN
                  CALL ORDENA1F2I(NTOPMAX,ATOP,ITOP,JTOP)
                END IF
              ELSE
                IF(B(J,I).GT.ATOP(1))THEN
                  ATOP(1)=B(J,I)
                  ITOP(1)=I
                  JTOP(1)=J
                  CALL ORDENA1F2I(NTOPMAX,ATOP,ITOP,JTOP)
                END IF
              END IF
            END DO
          END DO
          WRITE(*,*)
        END IF
C------------------------------------------------------------------------------
        BMINF=B(NC1,NS1)
        BMAXF=BMINF
        DO I=NS1,NS2
          DO J=NC1,NC2
            IF(B(J,I).LT.BMINF) BMINF=B(J,I)
            IF(B(J,I).GT.BMAXF) BMAXF=B(J,I)
          END DO
        END DO
20      WRITE(CDUMMY,*)BMINF
        WRITE(*,100)'Minimum value '
        BMIN=READF(CDUMMY)
        WRITE(CDUMMY,*)BMAXF
        WRITE(*,100)'Maximum value '
        BMAX=READF(CDUMMY)
        WRITE(CDUMMY,*)NBINMAX
        WRITE(*,100)'No. of bins '
        NBIN=READILIM(CDUMMY,1,NBINMAX)
C------------------------------------------------------------------------------
        DO K=1,NBIN
          NY(K)=0
        END DO
        DB=(BMAX-BMIN)/REAL(NBIN)
        DO I=NS1,NS2
          DO J=NC1,NC2
            K=INT((B(J,I)-BMIN)/DB)+1
            IF((K.GE.1).AND.(K.LE.NBIN)) NY(K)=NY(K)+1
          END DO
        END DO
C
        YMAX=0.
        YMIN=-0.3
        DO K=1,NBIN
          IF(NY(K).GT.0)THEN
            Y(K)=ALOG10(REAL(NY(K)))
            IF(Y(K).GT.YMAX) YMAX=Y(K)
          ELSE
            Y(K)=-1.
          END IF
        END DO
        DY=YMAX-YMIN
        YMAX=YMAX+DY/50.
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
ccc         CALL RPGERAS
            CALL BUTTSPR(0.10,0.95,0.10,0.70)
            CALL RPGENV(BMIN,BMAX,YMIN,YMAX,0,0)
            IF(NTYPE.EQ.3)THEN
              CALL PGSCI(3)
              CALL PGSLS(2)
              CALL PGMOVE(ATOP(1),YMIN)
              CALL PGDRAW(ATOP(1),YMAX)
              CALL PGSLS(1)
              CALL PGSCI(1)
            END IF
            CALL RPGENV(1.,REAL(NBIN),YMIN,YMAX,0,-2)
          ELSE
            CALL PGENV(BMIN,BMAX,YMIN,YMAX,0,-2)
ccc         CALL PGWINDOW(BMIN,BMAX,YMIN,YMAX)
            CALL PGBOX('BCNITS',0.0,0,'BCNITS',0.0,0)
            IF(NTYPE.EQ.3)THEN
              CALL PGSCI(3)
              CALL PGSLS(2)
              CALL PGMOVE(ATOP(1),YMIN)
              CALL PGDRAW(ATOP(1),YMAX)
              CALL PGSLS(1)
              CALL PGSCI(1)
            END IF
            CALL PGWINDOW(1.,REAL(NBIN),YMIN,YMAX)
          END IF
          CALL PGIDEN_RED
        END DO
        DO K=1,NBIN
          X(K)=REAL(K)
        END DO
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
          CALL PGBIN(NBIN,X,Y,.TRUE.)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          IF(NTYPE.EQ.1)THEN
            CALL PGLABEL('No. of counts','Log(N\dpixels\u)',NAMFIL)
          ELSE
            CALL PGLABEL('Deviations of the mean (in \gs units)',
     +       'Log(N\dpixels\u)',NAMFIL)
          END IF
        END DO
C
        IF(NTYPE.EQ.3)THEN
          WRITE(*,100)'NOTE: the dashed green-line indicates the '
          WRITE(*,101)'TOP100 level'
        END IF
        WRITE(*,100)'Replot (y/n) '
        CREPLOT(1:1)=READC('n','yn')
        IF(CREPLOT.EQ.'y') GOTO 20
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE SUBPLOT(J,I,LNEXT,PINTA)
C Esta subrutina dibuja la region alrededor del punto (J,I) utilizando
C escala de grises (funcion PGGRAY/PGIMAG de PGPLOT).
C En cada dibujo se representan NPP x NPP pixels.
C NOTA: NPP ha de ser impar.
C Un parametro importante es SIGTHRESHOLD, que determina el numero de veces
C sigma con que dibujamos el segundo plot de la region con C.R., y ademas
C es el valor utilizado para presenhalizar los pixeles candidatos a rayo
C cosmico.
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER I,J,L
        LOGICAL LNEXT,PINTA
C
        INTEGER NPPMAX
        INTEGER NBOTONES
        INTEGER NDPMAX                              !grado maximo del polinomio
        INTEGER NAUXMAX                   !numero maximo de imagenes auxiliares
        INTEGER MAXNCOEFF    !numero maximo de coeficientes (ajuste superficie)
        PARAMETER(NPPMAX=101)
        PARAMETER(NBOTONES=18)
        PARAMETER(NDPMAX=20)
        PARAMETER(NAUXMAX=3)          !si se cambia hacerlo en todos los sitios
        PARAMETER(MAXNCOEFF=100)
C
        INTEGER II,JJ,LL,K
        INTEGER III,JJJ
        INTEGER MII,MJJ
        INTEGER NS0,NS1,NC0,NS2,NC1,NC2
        INTEGER NPP
        INTEGER NPP2
        INTEGER NS,NC,NWIN
        INTEGER IXC,IYC
        INTEGER NFIT
        INTEGER NDP
        INTEGER II1,II2,JJ1,JJ2
        INTEGER NSKYLINE
        INTEGER CAUTOMODE
        INTEGER NB
        INTEGER NTOTCR
        INTEGER NTYPEPLOT
        INTEGER BMODE2(NBOTONES)
        INTEGER ANGLE3D
        INTEGER NAUXFRAME,NAUX,NCTOT0,NCTOT(NAUXMAX)
        INTEGER P,Q
        INTEGER ORDER(MAXNCOEFF),IOK,IPAR
        INTEGER NCOEFF,GX,GY
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL MEAN,SIGMA
        REAL A(NCMAX,NSMAX)                           !imagen normal de trabajo
        REAL AA(NCMAX,NSMAX)             !imagen backup para recuperar con UNDO
        REAL ERR(NCMAX,NSMAX)                                !imagen de errores
        REAL EERR(NCMAX,NSMAX)           !imagen backup para recuperar con UNDO
        REAL AREPLACE(NPPMAX,NPPMAX)
        REAL TR(6),FG,BG
        REAL FG3D,BG3D
        REAL SIGCR,FACTCR,MINVCR,SIGCRAUX
        REAL XC,YC
        REAL M(MAXNCOEFF,MAXNCOEFF),N(MAXNCOEFF)
        REAL SCALEROW(MAXNCOEFF),XSOL(MAXNCOEFF)
        REAL CCX1,CCX2,CCY1,CCY2
        REAL PERCENT1,PERCENT2
        REAL SIGTHRESHOLD
        REAL RQSKY
        REAL XPOL(NPPMAX),YPOL(NPPMAX),COEFP(NDPMAX+1),CHISQR
        REAL AUXFRAME(NCMAX,NSMAX,NAUXMAX)
        REAL FIXED_FG,FIXED_BG
        REAL FFACTOR
        REAL DEV
ccc     REAL X1VPORT,X2VPORT,Y1VPORT,Y2VPORT
        CHARACTER*1 CH,CFIX_FGBG,CERR,CCROS,CEXTRAPOL
        CHARACTER*20 LABEL2(NBOTONES)
        CHARACTER*70 GLABEL
        CHARACTER*50 CDUMMY
        CHARACTER*75 AUXFILE(NAUXMAX)
        CHARACTER*80 NAMFIL
        LOGICAL MASK(NPPMAX,NPPMAX),MASKSKY(NPPMAX)
        LOGICAL LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        LOGICAL CRFOUND
        LOGICAL CLEANAUTO,FIRSTTIME
        LOGICAL LHISTFILE
        LOGICAL GOODFIT,GOODFIT1,GOODFIT2
        LOGICAL LREAD
        LOGICAL MM3D
        LOGICAL CHANYES(NCMAX),SCANYES(NSMAX)
        LOGICAL LBEXIST
        LOGICAL IFPIXEL(NPPMAX,NPPMAX)
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL TOP100CLEAN(NCMAX,NSMAX)
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLKERR/ERR
        COMMON/BLK1/MEAN,SIGMA
        COMMON/BLK2/NS,NC,NWIN
        COMMON/BLK3A/LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        COMMON/BLK3B/SIGCR,FACTCR,MINVCR,SIGCRAUX
        COMMON/BLK5/NFIT,NPP,NDP
        COMMON/BLK4/MASK,NC1,NC2,NS1,NS2
        COMMON/BLK6/PERCENT1,PERCENT2
        COMMON/BLK7A/SIGTHRESHOLD
        COMMON/BLK7B/RQSKY
        COMMON/BLK8/MASKSKY
        COMMON/BLK9/CLEANAUTO,CAUTOMODE,CCROS
        COMMON/BLK10/NTOTCR
        COMMON/BLK11/LABEL2,BMODE2
        COMMON/BLK12/NAMFIL
        COMMON/BLK13/LHISTFILE
        COMMON/BLK15/NTYPEPLOT
        COMMON/BLK17A/NAUXFRAME
        COMMON/BLK17B/AUXFRAME
        COMMON/BLK17C/AUXFILE
        COMMON/BLK17D/STWV,DISP
        COMMON/BLK19/CHANYES,SCANYES
        COMMON/BLK20A/CFIX_FGBG
        COMMON/BLK20B/FIXED_FG,FIXED_BG
        COMMON/BLK21/ANGLE3D
        COMMON/BLK22/CERR
        COMMON/BLK23/AREPLACE
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
        COMMON/BLK3TOP100/TOP100CLEAN
C------------------------------------------------------------------------------
C do not clean automatically C.R. very near to the image border
        IF(CLEANAUTO)THEN
          IF((J.LE.NFIT).OR.(J.GE.NC-NFIT).OR.(I.LE.NFIT).OR.
     +     (I.GE.NS-NFIT))THEN
             WRITE(*,'(A,I5,I5,A)')'C.R. at:',J,I,' DO NOT removed.'
             IF(LHISTFILE)
     +        WRITE(80,'(A,I5,I5,A)')'C.R. at:',J,I,' DO NOT removed.'
             RETURN
          END IF
        END IF
C
        LNEXT=.TRUE.
        CRFOUND=.FALSE.
        MM3D=.FALSE.                                           !modo automatico
        IF(.NOT.CLEANAUTO)THEN
          IF((NTYPEPLOT.EQ.1).OR.(NTYPEPLOT.EQ.2))THEN
            CALL BUTTON(13,LABEL2(13),0)
            CALL BUTTON(14,LABEL2(14),0)
            CALL BUTTON(15,LABEL2(15),0)
          END IF
          IF(NAUXFRAME.GT.0) CALL BUTTON(16,LABEL2(16),0)
        END IF
C initialize the transformation matrix between array grid and world
C coordinates
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
        FIRSTTIME=.TRUE.
C
        NPP2=(NPP-1)/2
        NS1=I-NPP2
        NS2=I+NPP2
        NC1=J-NPP2
        NC2=J+NPP2
        IF(NS1.LT.1)THEN
          NS1=1
          NS2=NPP
        END IF
        IF (NS2.GT.NS)THEN
          NS1=NS-NPP+1
          NS2=NS
        END IF
        IF (NC1.LT.1)THEN
          NC1=1
          NC2=NPP
        END IF
        IF (NC2.GT.NC)THEN
          NC1=NC-NPP+1
          NC2=NC
        END IF
C
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            AA(JJ,II)=A(JJ,II)
          END DO
        END DO
        IF(CERR.EQ.'y')THEN
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              EERR(JJ,II)=ERR(JJ,II)
            END DO
          END DO
        END IF
C
10      FG=MEAN+3.*SIGMA
        BG=MEAN-3.*SIGMA
        IF(FG.GT.FG3D) FG=FG3D
        IF(BG.LT.BG3D) BG=BG3D
        IF(PINTA)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              CALL BUTTSPR(0.10,0.50,0.10,0.80)
            ELSE
              CALL PGPAGE
              CALL PGVPORT(0.10,0.50,0.10,0.80)
            END IF
          END DO
          IF(NTYPEPLOT.EQ.1)THEN
            CALL PLOT3DBARS(NC1,NC2,NS1,NS2,0,ANGLE3D,MM3D,BG3D,FG3D)
          ELSEIF(NTYPEPLOT.EQ.2)THEN
            CALL PLOT3DBARS(NC1,NC2,NS1,NS2,1,ANGLE3D,MM3D,BG3D,FG3D)
          ELSEIF(NTYPEPLOT.EQ.3)THEN
            DO ITERM=NTERM,1,-1
              IF(ITERM.EQ.1)THEN
                CALL RPGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +           REAL(NS2)+.6,1,0)
              ELSE
                CALL PGWNAD(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +           REAL(NS2)+.6)
                CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
              END IF
              CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
              CALL PGMTEXT('T',1.,0.,0.,'file: '//NAMFIL)
            END DO
          END IF
        END IF
C
C calculamos un nuevo valor de MEAN y SIGMA, usando solo la fraccion de
C pixels con senhal indicada a traves de PERCENT1 y PERCENT2
        CALL REFINE(NS1,NS2,NC1,NC2)
C
        FG=MEAN+2.*SIGTHRESHOLD*SIGMA
        BG=MEAN-2.*SIGTHRESHOLD*SIGMA
        IF((NTYPEPLOT.EQ.1).OR.(NTYPEPLOT.EQ.2))THEN
          IF(FG.GT.FG3D) FG=FG3D
          IF(BG.LT.BG3D) BG=BG3D
        END IF
        IF(PINTA)THEN
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            IF(ITERM.EQ.1)THEN
              CALL BUTTSPR(0.55,0.95,0.10,0.80)
              CALL RPGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +         REAL(NS2)+.6,1,0)
            ELSE
              CALL PGVPORT(0.55,0.95,0.10,0.80)
              CALL PGWNAD(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +         REAL(NS2)+.6)
              CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
            END IF
            WRITE(GLABEL,'(A,I4)')'Channel: ',J
            CALL PGMTEXT('B',3.,0.,0.,GLABEL)
            WRITE(GLABEL,'(A,I4)')'Scan: ',I
            CALL PGMTEXT('B',3.,1.,1.,GLABEL)
            CALL PGIDEN_RED
          END DO
        END IF
        CONTINUE
C
C Imponemos que los pixeles que constituyen el C.R. formen un conjunto
C conexo (subroutine WALKER).
C MASKSKY() indica si hemos encontrado alguna linea de cielo en el canal JJ,
C por lo que ahora hay que inicializarlo a FALSE. Mas adelante se testea la
C presencia de lineas de cielo.
        DO JJ=1,NPP
          MASKSKY(JJ)=.FALSE.
        END DO
        CALL WALKER(J,I,.FALSE.)
C
C Testeamos si la busqueda automatica ha podido encontrar alguna linea
C de cielo (el RQSKY% de los pixels de una columna presenhalizados). En dicho
C caso, eliminamos la presenhalizacion de todos los puntos de la citada
C columna por si acaso (el usuario tiene que decidir si algun pixel
C sobre la linea de cielo es o no rayo cosmico). Si el rayo cosmico es
C ENORME y vertical, podria darse el caso de que el programa lo confundiese
C con una linea de cielo, aunque para eso esta el usuario. Si se detecta
C una posible linea de cielo (y se elimina) hay que volver a chequear que
C todos los puntos presenhalizados forman un conjunto conexo.
C MASKSKY() determina si en el canal JJ se ha encontrado linea de cielo
        DO JJ=NC1,NC2
          NSKYLINE=0
          DO II=NS1,NS2
            IF(MASK(JJ-NC1+1,II-NS1+1))NSKYLINE=NSKYLINE+1
          END DO
          IF(NSKYLINE.GE.INT(RQSKY*REAL(NPP)/100.))THEN
            WRITE(*,'(A,I4)')'WARNING: Possible Sky Line '//
     +       'detected at channel ',JJ
            MASKSKY(JJ-NC1+1)=.TRUE.
            CALL WALKER(J,I,.FALSE.)
c El siguiente codigo (en minusculas) evita que la eliminacion de una linea
c de cielo de lugar a que WALKER seleccione otra linea de cielo (situada mas
c a la izquierda).
            do jjj=nc1,nc2
              if(.not.masksky(jjj-nc1+1))then
                nskyline=0
                do iii=ns1,ns2
                  if(mask(jjj-nc1+1,iii-ns1+1)) nskyline=nskyline+1
                end do
                if(nskyline.ge.int(rqsky*real(npp)/100.))then
                  write(*,'(a,i4)')'warning: possible sky line '//
     +             'detected at channel ',jjj
                  masksky(jjj-nc1+1)=.true.
                  call walker(j,i,.false.)
                end if
              end if
            end do
c
          END IF
        END DO
C
        IF(PINTA)THEN
          IF(CFIX_FGBG.EQ.'y')THEN
            FG=FIXED_FG
            BG=FIXED_BG
          END IF
          IF(FG.EQ.BG)THEN
            IF(FG.EQ.0.)THEN
              FG=1.
              BG=-1.
            ELSE
              FG=1.1*FG
              BG=0.9*BG
            END IF
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          END DO
        END IF
        CALL WALKER(J,I,.TRUE.)
C El usuario ahora indica que pixels habran de ser interpolados, y de
C que modo.
C------------------------------------------------------------------------------
20      CONTINUE
        IF(CLEANAUTO)THEN
          IF(FIRSTTIME)THEN
            NB=-1                                                 !ningun boton
            IF(CAUTOMODE.EQ.1)THEN
              CH='x'
            ELSEIF(CAUTOMODE.EQ.2)THEN
              CH='y'
            ELSEIF(CAUTOMODE.EQ.3)THEN
              CH='a'
            ELSE
              CH='z'
            END IF
            FIRSTTIME=.FALSE.
          ELSE
            CH='c'
          END IF
        ELSE
          IF(FIRSTTIME)THEN
            FIRSTTIME=.FALSE.
          END IF
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
        END IF
C------------------------------------------------------------------------------
        IF(CH.EQ.'?')THEN
          CALL DISPLAYMENU
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.12).OR.(CH.EQ.'e'))THEN
          LNEXT=.FALSE.
          IF(LHISTFILE)THEN
            IF(CRFOUND) WRITE(80,'(A,I5,I5)')'C.R. at:',J,I
          END IF
          WRITE(*,*)
          WRITE(*,104)'--[INFORMATION]--> NO. OF CR. REMOVED: ',
     +     NTOTCR
          WRITE(*,100)' current image: '
          WRITE(*,101)NAMFIL(1:TRUELEN(NAMFIL))
          RETURN
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.5).OR.(CH.EQ.'s').OR.(CH.EQ.',').OR.(CH.EQ.'/'))THEN
          CALL BUTTON(5,LABEL2(5),5)
          IF(CH.EQ.',')THEN
            FG=FG3D
            BG=BG3D
          ELSEIF(CH.EQ.'/')THEN
            CALL ZSCALE(A,NC1,NC2,NS1,NS2,BG,FG)
          ELSE
            WRITE(*,101)'REMEMBER: <,> set limits to YMIN,YMAX'
            WRITE(*,100)'Background    : '
            WRITE(*,*)BG
            WRITE(*,100)'Foreground    : '
            WRITE(*,*)FG
            WRITE(*,100)'New background '
            WRITE(CDUMMY,*)BG
            BG=READF(CDUMMY)
            WRITE(*,100)'New foreground '
            WRITE(CDUMMY,*)FG
            FG=READF(CDUMMY)
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          END DO
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1)) CALL MARCA(JJ,II)
            END DO
          END DO
          CALL BUTTON(5,LABEL2(5),0)
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.11).OR.(CH.EQ.'c'))THEN
          IF(PINTA) CALL BUTTON(11,LABEL2(11),5)
          WRITE(*,*)
          WRITE(*,104)'--[INFORMATION]--> NO. OF CR. REMOVED: ',
     +     NTOTCR
          WRITE(*,100)' current image: '
          WRITE(*,101)NAMFIL(1:TRUELEN(NAMFIL))
          IF(PINTA)THEN
            WRITE(*,*)
            CALL BUTTON(11,LABEL2(11),0)
          END IF
          GOTO 30
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.17).OR.(CH.EQ.'k'))THEN
          CALL BUTTON(17,LABEL2(17),5)
          CHANYES(J)=.FALSE.
          CALL BUTTON(17,LABEL2(17),0)
          GOTO 30
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.18).OR.(CH.EQ.'i'))THEN
          CALL BUTTON(18,LABEL2(18),5)
          SCANYES(I)=.FALSE.
          CALL BUTTON(18,LABEL2(18),0)
          GOTO 30
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.7).OR.(CH.EQ.'r'))THEN
          CALL BUTTON(7,LABEL2(7),5)
          WRITE(*,101)'removing crosses...'
          DO II=1,NPP
            DO JJ=1,NPP
              MASK(JJ,II)=.FALSE.
            END DO
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          END DO
          WRITE(*,*)
          CALL BUTTON(7,LABEL2(7),0)
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.8).OR.(CH.EQ.'t'))THEN
          CALL BUTTON(8,LABEL2(8),5)
          WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL BUTTON(6,CDUMMY(1:L),3)
          DO II=1,NBOTONES
            IF(BMODE2(II).GE.0)THEN
              IF((II.NE.6).AND.(II.NE.8))THEN
                CALL BUTTON(II,LABEL2(II),3)
              END IF
            END IF
          END DO
214       WRITE(*,101)'Press mouse buttom...'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
          IF((NB.NE.8).AND.(CH.NE.'t'))THEN
            IF((CH.EQ.'A').OR.(CH.EQ.'D').OR.(CH.EQ.'X')) THEN
              IXC=INT(XC+0.5)
              IYC=INT(YC+0.5)
              IF((IXC.GE.NC1).AND.(IXC.LE.NC2).AND.(IYC.GE.NS1).AND.
     +         (IYC.LE.NS2))THEN
                WRITE(*,100)'Position   : '
                WRITE(*,*)IXC,IYC
                WRITE(*,100)'Pixel value: '
                WRITE(*,*)A(IXC,IYC)
                WRITE(*,100)'MEAN,SIGMA : '
                WRITE(*,*)MEAN,SIGMA
                WRITE(*,100)'SIGMA times over MEAN: '
                WRITE(*,*)(A(IXC,IYC)-MEAN)/SIGMA
                WRITE(*,*)
              ELSE
                WRITE(*,101)'ERROR: Cursor outside plot region.'
              END IF
            ELSE
              WRITE(*,101)'ERROR: No mouse buttom detected.'
            END IF
            GOTO 214
          END IF
          WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL BUTTON(6,CDUMMY(1:L),0)
          DO II=1,NBOTONES
            IF(II.NE.6)THEN
              IF(BMODE2(II).GE.0) CALL BUTTON(II,LABEL2(II),0)
            END IF
          END DO
          IF(.NOT.CRFOUND) CALL BUTTON(9,LABEL2(9),3)
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.9).OR.(CH.EQ.'u'))THEN
          CALL BUTTQEX(9,LBEXIST)
          IF(LBEXIST)THEN
            CALL BUTTON(9,LABEL2(9),5)
            WRITE(*,117)'Restoring Area around ',J,I
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                A(JJ,II)=AA(JJ,II)
                TOP100CLEAN(JJ,II)=.FALSE.
              END DO
            END DO
            IF(CERR.EQ.'y')THEN
              DO II=NS1,NS2
                DO JJ=NC1,NC2
                  ERR(JJ,II)=EERR(JJ,II)
                END DO
              END DO
            END IF
            NTOTCR=NTOTCR-1
            CRFOUND=.FALSE.
            CALL BUTTON(9,LABEL2(9),0)
            CALL BUTTON(9,LABEL2(9),3)
            GOTO 10
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.13).OR.(CH.EQ.'+'))THEN
          CALL BUTTQEX(13,LBEXIST)
          IF(LBEXIST)THEN
            CALL BUTTON(13,LABEL2(13),5)
            ANGLE3D=ANGLE3D+90
            IF(ANGLE3D.EQ.360) ANGLE3D=0
            IF(NTYPEPLOT.EQ.1)THEN
              CALL PLOT3DBARS(NC1,NC2,NS1,NS2,0,ANGLE3D,MM3D,BG3D,FG3D)
            ELSEIF(NTYPEPLOT.EQ.2)THEN
              CALL PLOT3DBARS(NC1,NC2,NS1,NS2,1,ANGLE3D,MM3D,BG3D,FG3D)
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(ITERM.EQ.1)THEN
                CALL BUTTSPR(0.55,0.95,0.10,0.80)
                CALL RPGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +           REAL(NS2)+.6,1,0)
              ELSE
                CALL PGVPORT(0.55,0.95,0.10,0.80)
                CALL PGWNAD(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +           REAL(NS2)+.6)
ccc                CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
              END IF
            END DO
            CALL BUTTON(13,LABEL2(13),0)
            GOTO 20
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.14).OR.(CH.EQ.'-'))THEN
          CALL BUTTQEX(14,LBEXIST)
          IF(LBEXIST)THEN
            CALL BUTTON(14,LABEL2(14),5)
            ANGLE3D=ANGLE3D-90
            IF(ANGLE3D.LT.0) ANGLE3D=ANGLE3D+360
            IF(NTYPEPLOT.EQ.1)THEN
              CALL PLOT3DBARS(NC1,NC2,NS1,NS2,0,ANGLE3D,MM3D,BG3D,FG3D)
            ELSEIF(NTYPEPLOT.EQ.2)THEN
              CALL PLOT3DBARS(NC1,NC2,NS1,NS2,1,ANGLE3D,MM3D,BG3D,FG3D)
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(ITERM.EQ.1)THEN
                CALL BUTTSPR(0.55,0.95,0.10,0.80)
                CALL RPGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +           REAL(NS2)+.6,1,0)
              ELSE
                CALL PGVPORT(0.55,0.95,0.10,0.80)
                CALL PGWNAD(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +           REAL(NS2)+.6)
ccc                CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
              END IF
            END DO
            CALL BUTTON(14,LABEL2(14),0)
            GOTO 20
          END IF
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.15)THEN
          IF(MM3D)THEN
            LABEL2(15)='3D(auto)'
            MM3D=.FALSE.
          ELSE
            LABEL2(15)='3D(fixed)'
            MM3D=.TRUE.
          END IF
          CALL BUTTON(15,LABEL2(15),0)
          GOTO 20
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.16).OR.(CH.EQ.'f'))THEN
          CALL BUTTQEX(16,LBEXIST)
          IF(LBEXIST)THEN
            CALL BUTTON(16,LABEL2(16),1)
            WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL BUTTON(6,CDUMMY(1:L),3)
            DO II=1,NBOTONES
              IF(BMODE2(II).GE.0)THEN
                IF((II.NE.6).AND.(II.NE.16))THEN
                  CALL BUTTON(II,LABEL2(II),3)
                END IF
              END IF
            END DO
            IF(CRFOUND)THEN
              WRITE(*,101)'WARNING: this action will '//
     +         'restore cleaned region'
              WRITE(*,100)'Press <CR> to continue...'
              READ(*,*)
            END IF
            NCTOT0=0
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                NCTOT0=NCTOT0+AA(JJ,II)
              END DO
            END DO
            DO NAUX=1,NAUXFRAME
              NCTOT(NAUX)=0
              DO II=NS1,NS2
                DO JJ=NC1,NC2
                  NCTOT(NAUX)=NCTOT(NAUX)+AUXFRAME(JJ,II,NAUX)
                END DO
              END DO
            END DO
            WRITE(*,*)
            WRITE(*,101)'* Press any key to change frame'
            WRITE(*,101)'* Press <q>, <x> or <right mouse '//
     +       'button> to exit'
            WRITE(*,*)
            CALL CHECK_AUXFRAME(J,I,DEV)
            NAUX=NAUXFRAME
212         NAUX=NAUX+1
            IF(NAUX.EQ.NAUXFRAME+1)THEN
              WRITE(*,101)'* Displayed image is Original Frame'
              WRITE(*,101)'- '//NAMFIL(1:TRUELEN(NAMFIL))
              DO II=NS1,NS2
                DO JJ=NC1,NC2
                  A(JJ,II)=AA(JJ,II)
                END DO
              END DO
            ELSEIF(NAUX.EQ.NAUXFRAME+2)THEN
              WRITE(*,101)'* Displayed image is Mean of Aux frames'
              DO II=NS1,NS2
                DO JJ=NC1,NC2
                  A(JJ,II)=AREPLACE(JJ-NC1+1,II-NS1+1)
                END DO
              END DO
              NAUX=0
            ELSE
              WRITE(*,110)'* Displayed image is auxiliary frame #',NAUX
              WRITE(*,101)'- '//AUXFILE(NAUX)
              DO II=NS1,NS2
                DO JJ=NC1,NC2
ccc                  A(JJ,II)=AUXFRAME(JJ,II,NAUX)*NCTOT0/NCTOT(NAUX)
                  A(JJ,II)=AUXFRAME(JJ,II,NAUX)
                END DO
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
            END DO
            CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
            CALL CHLOWER(CH)
            IF((CH.NE.'q').AND.(CH.NE.'x')) GOTO 212
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                A(JJ,II)=AA(JJ,II)
              END DO
            END DO
            WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            CALL BUTTON(6,CDUMMY(1:L),0)
            DO II=1,NBOTONES
              IF(II.NE.6)THEN
                IF(BMODE2(II).GE.0) CALL BUTTON(II,LABEL2(II),0)
              END IF
            END DO
            IF(CRFOUND)THEN
              NTOTCR=NTOTCR-1
              CALL BUTTON(9,LABEL2(9),3)
              CRFOUND=.FALSE.
            END IF
            WRITE(*,101)'* Displayed image is Original Frame'
            WRITE(*,101)'- '//NAMFIL(1:TRUELEN(NAMFIL))
            GOTO 10
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.6).OR.(CH.EQ.'n'))THEN
          WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL BUTTON(6,CDUMMY(1:L),5)
          LREAD=.TRUE.
          DO WHILE((NFIT.LT.1).OR.LREAD)
            LREAD=.FALSE.
            WRITE(*,100)'No. of pixels for fit '
            WRITE(CDUMMY,*)NFIT
            NFIT=READI(CDUMMY)
            IF(NFIT.LT.1)WRITE(*,101)'ERROR: invalid number. Try again.'
          END DO
          LREAD=.TRUE.
          DO WHILE(((NDP.LT.1).OR.(NDP.GT.NDPMAX)).OR.LREAD)
            LREAD=.FALSE.
            WRITE(*,100)'Polynomial degree for fit '
            WRITE(CDUMMY,*)NDP-1
            NDP=READILIM(CDUMMY,0,2*NFIT-1)
            NDP=NDP+1
            IF((NDP.LT.1).OR.(NDP.GT.NDPMAX))
     +       WRITE(*,101)'ERROR: invalid number. Try again.'
          END DO
          WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL BUTTON(6,CDUMMY(1:L),0)
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.1).OR.(CH.EQ.'x'))THEN
          IF(PINTA) CALL BUTTON(1,LABEL2(1),5)
          GOODFIT=.FALSE.                  !comprobamos que hay cruces marcadas
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT=.TRUE.
            END DO
          END DO
          IF(.NOT.GOODFIT)THEN
            WRITE(*,101)'ERROR: no. of pixels selected=0.'
            IF(PINTA) CALL BUTTON(1,LABEL2(1),0)
            GOTO 20
          END IF
          WRITE(*,117)'Cleaning C.R. around pixel ',J,I
          GOODFIT1=.TRUE.!chequeamos que el numero de puntos a cada lado existe
          GOODFIT2=.TRUE.
          DO II=NS1,NS2
            DO JJ=NC1,NC1+NFIT-1
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT1=.FALSE.
            END DO
            DO JJ=NC2-NFIT+1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT2=.FALSE.
            END DO
          END DO
          GOODFIT=(GOODFIT1.AND.GOODFIT2)
          IF(.NOT.GOODFIT)THEN
            IF((GOODFIT1.OR.GOODFIT2).AND.(.NOT.CLEANAUTO))THEN
              WRITE(*,101)'WARNING: some pixel(s) is(are) too near to'//
     +         ' the window edge(s).'
              WRITE(*,100)'Do you want to extrapolate the fit (y/n) '
              CEXTRAPOL(1:1)=READC('y','yn')
              IF(CEXTRAPOL.EQ.'n')THEN
                IF(PINTA) CALL BUTTON(1,LABEL2(1),0)
                GOTO 20
              ELSE
                IF(NFIT.LT.NDP)THEN
                  WRITE(*,101)'ERROR: no. of pixels is less than '//
     +             'polynomial degree + 1'
                  IF(PINTA) CALL BUTTON(1,LABEL2(1),0)
                  GOTO 20
                END IF
              END IF
            ELSE
              WRITE(*,101)'ERROR: some pixel(s) is(are) too near to '//
     +         'the window edge(s).'
              IF(PINTA) CALL BUTTON(1,LABEL2(1),0)
              GOTO 20
            END IF
          END IF
          CRFOUND=.TRUE.
          DO II=NS1,NS2
            JJ=NC1                      !buscamos primer canal por la izquierda
            DO WHILE((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(JJ.LT.NC2))
              JJ=JJ+1
            END DO
            IF((JJ.LT.NC2).OR.(MASK(NC2-NC1+1,II-NS1+1)))THEN
              JJ1=JJ-1             !JJ es el primer canal (del scan II) marcado
              JJ=NC2                      !buscamos ultimo canal por la derecha
              DO WHILE(.NOT.MASK(JJ-NC1+1,II-NS1+1).AND.(JJ.GT.NC1))
                JJ=JJ-1
              END DO
              JJ2=JJ+1             !JJ es el ultimo canal (del scan II) marcado
              LL=0       !metemos en XPOL(),YPOL() los puntos para el polinomio
              DO MJJ=JJ1-NFIT+1,JJ1
                IF(MJJ.GE.NC1)THEN
                  LL=LL+1
                  XPOL(LL)=REAL(MJJ-NC1+1)
                  YPOL(LL)=A(MJJ,II)
                END IF
              END DO
              DO MJJ=JJ2,JJ2+NFIT-1
                IF(MJJ.LE.NC2)THEN
                  LL=LL+1
                  XPOL(LL)=REAL(MJJ-NC1+1)
                  YPOL(LL)=A(MJJ,II)
                END IF
              END DO
              CALL POLFIT(XPOL,YPOL,YPOL,LL,NDP,0,COEFP,CHISQR)
              DO MJJ=JJ1+1,JJ2-1             !sustituimos valores del polinomio
                A(MJJ,II)=COEFP(NDP)                   !calculamos el polinomio
                DO LL=NDP-1,1,-1
                  A(MJJ,II)=A(MJJ,II)*REAL(MJJ-NC1+1)+COEFP(LL)
                END DO
                MASK(MJJ-NC1+1,II-NS1+1)=.FALSE.
                TOP100CLEAN(MJJ,II)=.TRUE.
              END DO
              IF(CERR.EQ.'y')THEN
                LL=0     !metemos en XPOL(),YPOL() los puntos para el polinomio
                DO MJJ=JJ1-NFIT+1,JJ1
                  IF(MJJ.GE.NC1)THEN
                    LL=LL+1
                    XPOL(LL)=REAL(MJJ-NC1+1)
                    YPOL(LL)=ERR(MJJ,II)
                  END IF
                END DO
                DO MJJ=JJ2,JJ2+NFIT-1
                  IF(MJJ.LE.NC2)THEN
                    LL=LL+1
                    XPOL(LL)=REAL(MJJ-NC1+1)
                    YPOL(LL)=ERR(MJJ,II)
                  END IF
                END DO
                CALL POLFIT(XPOL,YPOL,YPOL,LL,NDP,0,COEFP,CHISQR)
                DO MJJ=JJ1+1,JJ2-1           !sustituimos valores del polinomio
                  ERR(MJJ,II)=COEFP(NDP)               !calculamos el polinomio
                  DO LL=NDP-1,1,-1
                    ERR(MJJ,II)=ERR(MJJ,II)*REAL(MJJ-NC1+1)+COEFP(LL)
                  END DO
                END DO
              END IF
            END IF
          END DO
          NTOTCR=NTOTCR+1
          IF(PINTA)THEN
            CALL BUTTON(9,LABEL2(9),0)
            CALL BUTTON(1,LABEL2(1),0)
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.2).OR.(CH.EQ.'y'))THEN
          IF(PINTA) CALL BUTTON(2,LABEL2(2),5)
          GOODFIT=.FALSE.                  !comprobamos que hay cruces marcadas
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT=.TRUE.
            END DO
          END DO
          IF(.NOT.GOODFIT)THEN
            WRITE(*,101)'ERROR: no. of pixels selected=0.'
            IF(PINTA) CALL BUTTON(2,LABEL2(2),0)
            GOTO 20
          END IF
          WRITE(*,117)'Cleaning C.R. around pixel ',J,I
          GOODFIT1=.TRUE.!chequeamos que el numero de puntos a cada lado existe
          GOODFIT2=.TRUE.
          DO JJ=NC1,NC2
            DO II=NS1,NS1+NFIT-1
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT1=.FALSE.
            END DO
            DO II=NS2-NFIT+1,NS2
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT2=.FALSE.
            END DO
          END DO
          GOODFIT=(GOODFIT1.AND.GOODFIT2)
          IF(.NOT.GOODFIT)THEN
            IF((GOODFIT1.OR.GOODFIT2).AND.(.NOT.CLEANAUTO))THEN
              WRITE(*,101)'WARNING: some pixel(s) is(are) too near to'//
     +         ' the window edge(s).'
              WRITE(*,100)'Do you want to extrapolate the fit (y/n) '
              CEXTRAPOL(1:1)=READC('y','yn')
              IF(CEXTRAPOL.EQ.'n')THEN
                IF(PINTA) CALL BUTTON(2,LABEL2(2),0)
                GOTO 20
              ELSE
                IF(NFIT.LT.NDP)THEN
                  WRITE(*,101)'ERROR: no. of pixels is less than '//
     +             'polynomial degree + 1'
                  IF(PINTA) CALL BUTTON(2,LABEL2(2),0)
                  GOTO 20
                END IF
              END IF
            ELSE
              WRITE(*,101)'ERROR: some pixel(s) is(are) too near to '//
     +         'the window edge(s).'
              IF(PINTA) CALL BUTTON(2,LABEL2(2),0)
              GOTO 20
            END IF
          END IF
          CRFOUND=.TRUE.
          DO JJ=NC1,NC2
            II=NS1                                        !buscamos primer scan
            DO WHILE((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(II.LT.NS2))
              II=II+1
            END DO
            IF((II.LT.NS2).OR.(MASK(JJ-NC1+1,NS2-NS1+1)))THEN
              II1=II-1            !II1 es el primer scan (del canal JJ) marcado
              II=NS2                                      !buscamos ultimo scan
              DO WHILE(.NOT.MASK(JJ-NC1+1,II-NS1+1).AND.(II.GT.NS1))
                II=II-1
              END DO
              II2=II+1            !II2 es el ultimo scan (del canal JJ) marcado
              LL=0       !metemos en XPOL(),YPOL() los puntos para el polinomio
              DO MII=II1-NFIT+1,II1
                IF(MII.GE.NS1)THEN
                  LL=LL+1
                  XPOL(LL)=REAL(MII-NS1+1)
                  YPOL(LL)=A(JJ,MII)
                END IF
              END DO
              DO MII=II2,II2+NFIT-1
                IF(MII.LE.NS2)THEN
                  LL=LL+1
                  XPOL(LL)=REAL(MII-NS1+1)
                  YPOL(LL)=A(JJ,MII)
                END IF
              END DO
              CALL POLFIT(XPOL,YPOL,YPOL,LL,NDP,0,COEFP,CHISQR)
              DO MII=II1+1,II2-1             !sustituimos valores del polinomio
                A(JJ,MII)=COEFP(NDP)                   !calculamos el polinomio
                DO LL=NDP-1,1,-1
                  A(JJ,MII)=A(JJ,MII)*REAL(MII-NS1+1)+COEFP(LL)
                END DO
                MASK(JJ-NC1+1,MII-NS1+1)=.FALSE.
                TOP100CLEAN(JJ,MII)=.TRUE.
              END DO
              IF(CERR.EQ.'y')THEN
                LL=0     !metemos en XPOL(),YPOL() los puntos para el polinomio
                DO MII=II1-NFIT+1,II1
                  IF(MII.GE.NS1)THEN
                    LL=LL+1
                    XPOL(LL)=REAL(MII-NS1+1)
                    YPOL(LL)=ERR(JJ,MII)
                  END IF
                END DO
                DO MII=II2,II2+NFIT-1
                  IF(MII.LE.NS2)THEN
                    LL=LL+1
                    XPOL(LL)=REAL(MII-NS1+1)
                    YPOL(LL)=ERR(JJ,MII)
                  END IF
                END DO
                CALL POLFIT(XPOL,YPOL,YPOL,LL,NDP,0,COEFP,CHISQR)
                DO MII=II1+1,II2-1           !sustituimos valores del polinomio
                  ERR(JJ,MII)=COEFP(NDP)               !calculamos el polinomio
                  DO LL=NDP-1,1,-1
                    ERR(JJ,MII)=ERR(JJ,MII)*REAL(MII-NS1+1)+COEFP(LL)
                  END DO
                END DO
              END IF
            END IF
          END DO
          NTOTCR=NTOTCR+1
          IF(PINTA)THEN
            CALL BUTTON(9,LABEL2(9),0)
            CALL BUTTON(2,LABEL2(2),0)
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.3).OR.(CH.EQ.'a'))THEN
          IF(PINTA) CALL BUTTON(3,LABEL2(3),5)
          GOODFIT=.FALSE.                  !comprobamos que hay cruces marcadas
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT=.TRUE.
            END DO
          END DO
          IF(.NOT.GOODFIT)THEN
            WRITE(*,101)'ERROR: no. of pixels selected=0.'
            IF(PINTA) CALL BUTTON(3,LABEL2(3),0)
            GOTO 20
          END IF
C
          WRITE(*,117)'Cleaning C.R. around pixel ',J,I
          CRFOUND=.TRUE.
C
          DO II=NS1,NS2               !definimos los puntos a usar en el ajuste
            DO JJ=NC1,NC2
              IFPIXEL(JJ-NC1+1,II-NS1+1)=.FALSE.
            END DO
          END DO
C
          DO II=NS1,NS2                                          !buscamos en X
            JJ=NC1                                 !comenzamos por la izquierda
            DO WHILE((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(JJ.LT.NC2))
              JJ=JJ+1
            END DO
            IF(JJ.LT.NC2)THEN
              JJ1=JJ-1
              DO JJ=JJ1-NFIT+1,JJ1
                IF(JJ.GE.NC1)THEN
                  IFPIXEL(JJ-NC1+1,II-NS1+1)=.TRUE.
                END IF
              END DO
              JJ=NC2                                 !comenzamos por la derecha
              DO WHILE((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(JJ.GT.NC1))
                JJ=JJ-1
              END DO
              JJ2=JJ+1
              DO JJ=JJ2,JJ2+NFIT-1
                IF(JJ.LE.NC2)THEN
                  IFPIXEL(JJ-NC1+1,II-NS1+1)=.TRUE.
                END IF
              END DO
            END IF
          END DO
C
          DO JJ=NC1,NC2                                          !buscamos en Y
            II=NS1                                        !comenzamos por abajo
            DO WHILE((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(II.LT.NS2))
              II=II+1
            END DO
            IF(II.LT.NS2)THEN
              II1=II-1
              DO II=II1-NFIT+1,II1
                IF(II.GE.NS1)THEN
                  IFPIXEL(JJ-NC1+1,II-NS1+1)=.TRUE.
                END IF
              END DO
              II=NS2                                     !comenzamos por arriba
              DO WHILE((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(II.GT.NS1))
                II=II-1
              END DO
              II2=II+1
              DO II=II2,II2+NFIT-1
                IF(II.LE.NS2)THEN
                  IFPIXEL(JJ-NC1+1,II-NS1+1)=.TRUE.
                END IF
              END DO
            END IF
          END DO
C normalizamos el recorrido de las variables X e Y al intervalo [-1,1]
          CCX1=2./REAL(NC2-NC1)
          CCX2=REAL(NC1+NC2)/REAL(NC2-NC1)
          DO JJ=NC1,NC2
            XPOL(JJ-NC1+1)=CCX1*REAL(JJ)-CCX2
          END DO
          CCY1=2./REAL(NS2-NS1)
          CCY2=REAL(NS1+NS2)/REAL(NS2-NS1)
          DO II=NS1,NS2
            YPOL(II-NS1+1)=CCY1*REAL(II)-CCY2
          END DO
C resolvemos un sistema del tipo M(II,JJ) * X = N(II)
          GX=NDP-1   !tomamos como grado en cada eje el de los ajustes en X e Y
          GY=NDP-1
          NCOEFF=(GX+1)*(GY+1)
          DO WHILE(NCOEFF.GT.MAXNCOEFF)
            GX=GX-1
            GY=GY-1
            NCOEFF=(GX+1)*(GY+1)
            WRITE(*,110)'WARNING: new polynomial degree in X and Y: ',GX
          END DO
          II=0
          DO P=0,GX
            DO Q=0,GY
              II=II+1
              JJ=0
              DO III=0,GX
                DO JJJ=0,GY
                  JJ=JJ+1
                  M(II,JJ)=0.
                  DO NC0=NC1,NC2
                    DO NS0=NS1,NS2
                      IF(IFPIXEL(NC0-NC1+1,NS0-NS1+1))THEN
                        FFACTOR=1.
                        IF(III+P.NE.0)
     +                   FFACTOR=FFACTOR*(XPOL(NC0-NC1+1)**(III+P))
                        IF(JJJ+Q.NE.0)
     +                   FFACTOR=FFACTOR*(YPOL(NS0-NS1+1)**(JJJ+Q))
                        M(II,JJ)=M(II,JJ)+FFACTOR
CCC hemos cambiado las 2 lineas de abajo por las 6 lineas de arriba para
CCC evitar la indeterminacion 0.**0, que debe ser 1 para que funcione bien
ccc                     M(II,JJ)=M(II,JJ)+(XPOL(NC0-NC1+1)**(III+P))*
ccc     +                                    (YPOL(NS0-NS1+1)**(JJJ+Q))
                      END IF
                    END DO
                  END DO
                END DO
              END DO
              N(II)=0.
              DO NC0=NC1,NC2
                DO NS0=NS1,NS2
                  IF(IFPIXEL(NC0-NC1+1,NS0-NS1+1))THEN
                    FFACTOR=A(NC0,NS0)
                    IF(P.NE.0)FFACTOR=FFACTOR*(XPOL(NC0-NC1+1)**(P))
                    IF(Q.NE.0)FFACTOR=FFACTOR*(YPOL(NS0-NS1+1)**(Q))
                    N(II)=N(II)+FFACTOR
CCC hemos cambiado las 2 lineas de abajo por las 4 lineas de arriba para
CCC evitar la indeterminacion 0.**0, que debe ser 1 para que funcione bien
ccc                 N(II)=N(II)+A(NC0,NS0)*(XPOL(NC0-NC1+1)**(P))*
ccc     +                                     (YPOL(NS0-NS1+1)**(Q))
                  END IF
                END DO
              END DO
            END DO
          END DO
C resolvemos el sistema de ecuaciones lineales mediante descomposicion LU
          CALL LUDCMP(M,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,IOK,IPAR)
          IF(IOK.GT.0)THEN
            IF(N(ORDER(IOK)).EQ.0.)THEN
              WRITE(*,101)'ERROR: the system is redundant in LUDCMP.'
            ELSE
              WRITE(*,101)'ERROR: the system is inconsistent in LUDCMP.'
            END IF
            IF(PINTA) CALL BUTTON(3,LABEL2(3),0)
            GOTO 20
          ELSEIF(IOK.LT.0)THEN
            IF(PINTA) CALL BUTTON(3,LABEL2(3),0)
            GOTO 20
          END IF
          CALL LUSOLV(M,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,N,XSOL)
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1))THEN
                A(JJ,II)=0.
                K=0
                DO III=0,GX
                  DO JJJ=0,GY
                    K=K+1
                    FFACTOR=XSOL(K)
                    IF(III.NE.0)
     +               FFACTOR=FFACTOR*(XPOL(JJ-NC1+1)**(III))
                    IF(JJJ.NE.0)
     +               FFACTOR=FFACTOR*(YPOL(II-NS1+1)**(JJJ))
                    A(JJ,II)=A(JJ,II)+FFACTOR
CCC hemos cambiado las 2 lineas de abajo por las 6 lineas de arriba para
CCC evitar la indeterminacion 0.**0, que debe ser 1 para que funcione bien
ccc                 A(JJ,II)=A(JJ,II)+XSOL(K)*(XPOL(JJ-NC1+1)**(III))*
ccc     +                                        (YPOL(II-NS1+1)**(JJJ))
                  END DO
                END DO
                TOP100CLEAN(JJ,II)=.TRUE.
              END IF
            END DO
          END DO
C errores
          IF(CERR.EQ.'y')THEN
            II=0
            DO P=0,GX
              DO Q=0,GY
                II=II+1
                JJ=0
                N(II)=0.
                DO NC0=NC1,NC2
                  DO NS0=NS1,NS2
                    IF(IFPIXEL(NC0-NC1+1,NS0-NS1+1))THEN
                      FFACTOR=ERR(NC0,NS0)
                      IF(P.NE.0)FFACTOR=FFACTOR*(XPOL(NC0-NC1+1)**(P))
                      IF(Q.NE.0)FFACTOR=FFACTOR*(YPOL(NS0-NS1+1)**(Q))
                      N(II)=N(II)+FFACTOR
CCC hemos cambiado las 2 lineas de abajo por las 4 lineas de arriba para
CCC evitar la indeterminacion 0.**0, que debe ser 1 para que funcione bien
ccc                   N(II)=N(II)+ERR(NC0,NS0)*(XPOL(NC0-NC1+1)**(P))*
ccc     +                                         (YPOL(NS0-NS1+1)**(Q))
                    END IF
                  END DO
                END DO
              END DO
            END DO
C resolvemos el sistema de ecuaciones (no hace falta repetir CALL LUDCMP)
            CALL LUSOLV(M,NCOEFF,MAXNCOEFF,ORDER,SCALEROW,N,XSOL)
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                IF(MASK(JJ-NC1+1,II-NS1+1))THEN
                  ERR(JJ,II)=0.
                  K=0
                  DO III=0,GX
                    DO JJJ=0,GY
                      K=K+1
                      FFACTOR=XSOL(K)
                      IF(III.NE.0)
     +                 FFACTOR=FFACTOR*(XPOL(JJ-NC1+1)**(III))
                      IF(JJJ.NE.0)
     +                 FFACTOR=FFACTOR*(YPOL(II-NS1+1)**(JJJ))
                      ERR(JJ,II)=ERR(JJ,II)+FFACTOR
CCC hemos cambiado las 3 lineas de abajo por las 6 lineas de arriba para
CCC evitar la indeterminacion 0.**0, que debe ser 1 para que funcione bien
ccc                   ERR(JJ,II)=ERR(JJ,II)+
ccc     +                 XSOL(K)*(XPOL(JJ-NC1+1)**(III))*
ccc     +                         (YPOL(II-NS1+1)**(JJJ))
                    END DO
                  END DO
                END IF
              END DO
            END DO
          END IF
C
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              MASK(JJ-NC1+1,II-NS1+1)=.FALSE.
            END DO
          END DO
          NTOTCR=NTOTCR+1
          IF(PINTA)THEN
            CALL BUTTON(9,LABEL2(9),0)
            CALL BUTTON(3,LABEL2(3),0)
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.4).OR.(CH.EQ.'z'))THEN
          IF(PINTA) CALL BUTTON(4,LABEL2(4),5)
          IF(NAUXFRAME.LT.1)THEN
            WRITE(*,101)'ERROR: no. of auxiliary images=0.'
            IF(PINTA) CALL BUTTON(4,LABEL2(4),0)
            GOTO 20
          END IF
          GOODFIT=.FALSE.                  !comprobamos que hay cruces marcadas
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1)) GOODFIT=.TRUE.
            END DO
          END DO
          IF(.NOT.GOODFIT)THEN
            WRITE(*,101)'ERROR: no. of pixels selected=0.'
            IF(PINTA) CALL BUTTON(4,LABEL2(4),0)
            GOTO 20
          END IF
          WRITE(*,117)'Cleaning C.R. around pixel ',J,I
          CALL CHECK_AUXFRAME(J,I,DEV) !calculamos imagen promedio de AUXFRAME
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              IF(MASK(JJ-NC1+1,II-NS1+1))THEN
                A(JJ,II)=AREPLACE(JJ-NC1+1,II-NS1+1)
              END IF
            END DO
          END DO
C
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              MASK(JJ-NC1+1,II-NS1+1)=.FALSE.
            END DO
          END DO
          NTOTCR=NTOTCR+1
          IF(PINTA)THEN
            CALL BUTTON(9,LABEL2(9),0)
            CALL BUTTON(4,LABEL2(4),0)
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.1).OR.(NB.EQ.2).OR.(NB.EQ.3).OR.(NB.EQ.4).OR.
     +   (CH.EQ.'x').OR.(CH.EQ.'y').OR.(CH.EQ.'a').OR.(CH.EQ.'z'))THEN
          IF(PINTA)THEN
            IF((NTYPEPLOT.EQ.1).OR.(NTYPEPLOT.EQ.2))THEN
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(ITERM.EQ.1)THEN
                  CALL BUTTSPR(0.10,0.50,0.10,0.80)
                  CALL PGVPORT(0.10,0.50,0.10,0.80)
                  CALL PGWINDOW(0.,1.,0.,1.)
                  CALL PGSFS(1)
                  CALL PGSCI(0)
                  CALL PGRECT(0.,1.,0.,1.)
                  CALL PGSCI(1)
                ELSE
                  CALL PGPAGE
                END IF
              END DO
              IF(NTYPEPLOT.EQ.1)THEN
                CALL PLOT3DBARS(NC1,NC2,NS1,NS2,0,ANGLE3D,MM3D,
     +           BG3D,FG3D)
              ELSEIF(NTYPEPLOT.EQ.2)THEN
                CALL PLOT3DBARS(NC1,NC2,NS1,NS2,1,ANGLE3D,MM3D,
     +           BG3D,FG3D)
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF(ITERM.EQ.1)THEN
                  CALL BUTTSPR(0.55,0.95,0.10,0.80)
                  CALL RPGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +             REAL(NS2)+.6,1,0)
                ELSE
                  CALL PGVPORT(0.55,0.95,0.10,0.80)
                  CALL PGWNAD(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +             REAL(NS2)+.6)
                  CALL PGBOX('BCNTSI',0.0,0,'BCNTSI',0.0,0)
                END IF
                CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
              END DO
            END IF
          END IF
        END IF
C------------------------------------------------------------------------------
        IF((NB.EQ.10).OR.(CH.EQ.'#'))THEN
          CALL BUTTON(10,LABEL2(10),5)
203       WRITE(*,101)'% of Pixels with the lowest signal '//
     +     '(in the plot area) to be ignored in the'
          WRITE(*,101)'calculation of more refined values '//
     +     'for MEAN and SIGMA.'
          WRITE(*,102)'(current=',PERCENT1,')'
          WRITE(*,100)'New value '
          WRITE(CDUMMY,*)PERCENT1
          PERCENT1=READF(CDUMMY)
          IF((PERCENT1.LE.0.).OR.(PERCENT1.GT.100.))THEN
            WRITE(*,101)'ERROR: Invalid number. Try again.'
            GOTO 203
          END IF
          WRITE(*,101)'----------------------'
C
204       WRITE(*,101)'% of Pixels with the highest signal '//
     +     '(in the plot area) to be ignored in the'
          WRITE(*,101)'calculation of more refined values '//
     +     'for MEAN and SIGMA.'
          WRITE(*,102)'(current=',PERCENT2,')'
          WRITE(*,100)'New value '
          WRITE(CDUMMY,*)PERCENT2
          PERCENT2=READF(CDUMMY)
          IF((PERCENT2.LE.PERCENT1).OR.(PERCENT2.GT.100.))THEN
            WRITE(*,101)'ERROR: Invalid number. Try again.'
            GOTO 204
          END IF
          WRITE(*,101)'----------------------'
C
205       WRITE(*,101)'SIGTHRESHOLD is the number of times that '//
     +     'a pixel value exceeds MEAN '
          WRITE(*,101)'in SIGMA units'
          WRITE(*,102)'(current=',SIGTHRESHOLD,')'
          WRITE(*,100)'New value '
          WRITE(CDUMMY,*)SIGTHRESHOLD
          SIGTHRESHOLD=READF(CDUMMY)
          IF(SIGTHRESHOLD.LE.0.)THEN
            WRITE(*,101)'ERROR: Invalid number. Try again.'
            GOTO 205
          END IF
C
          WRITE(*,101)'----------------------'
C
206       WRITE(*,101)'% of the number of pixels in a column that '//
     +     'must be misinterpreted as C.R.'
          WRITE(*,101)'to suspect the presence of a Sky Line'
          WRITE(*,102)'(current=',RQSKY,')'
          WRITE(*,100)'New value '
          WRITE(CDUMMY,*)RQSKY
          RQSKY=READF(CDUMMY)
          IF((RQSKY.LE.0.).OR.(RQSKY.GT.100.))THEN
            WRITE(*,101)'ERROR: Invalid number. Try again.'
            GOTO 206
          END IF
C
          CALL BUTTON(10,LABEL2(10),0)
        END IF
C------------------------------------------------------------------------------
C seleccion de pixel con el cursor (boton izquierdo del raton)
        IF(NB.NE.0) GOTO 20
        IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) GOTO 20
        IXC=INT(XC+0.5)
        IYC=INT(YC+0.5)
        IF((IXC.GE.NC1).AND.(IXC.LE.NC2))THEN
          IF((IYC.GE.NS1).AND.(IYC.LE.NS2))THEN
            IF(.NOT.MASK(IXC-NC1+1,IYC-NS1+1))THEN
              CALL MARCA(IXC,IYC)
              WRITE(*,117)'Cursor at ',IXC,IYC
              MASK(IXC-NC1+1,IYC-NS1+1)=.TRUE.
            ELSE
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGIMAG(A,NCMAX,NSMAX,IXC,IXC,IYC,IYC,FG,BG,TR)
              END DO
              MASK(IXC-NC1+1,IYC-NS1+1)=.FALSE.
            END IF
          END IF
        END IF
        GOTO 20
30      CONTINUE        
        IF(LHISTFILE)THEN
          IF(CRFOUND) WRITE(80,'(A,I5,I5)')'C.R. at:',J,I
        END IF
        RETURN
C
100     FORMAT(A,$)
101     FORMAT(A,$)
102     FORMAT(A,F7.3,A)
104     FORMAT(A,I4)
110     FORMAT(A,I6)
117     FORMAT(A,I5,2X,I5)
C
        END
C
C******************************************************************************
C
        SUBROUTINE REFINE(NS1,NS2,NC1,NC2)
C Calcula un nuevo valor de MEAN y SIGMA utilizando solo el una fraccion
C de pixels comprendidos entre PERCENT1 y PERCENT2 (ordenados por orden
C creciente de senhal) en el cuadrado de NPP x NPP pixels.
C LREMOVEX=.TRUE.: elimina forma en direccion X antes de calcular MEAN,SIGMA
C LREMOVEY=.TRUE.: elimina forma en direccion Y antes de calcular MEAN,SIGMA
C En esta subrutina no se tiene en cuenta si se ha eliminado o no la forma
C promedio en X e Y antes de buscar el rayo cosmico (si se introduce dicha
C opcion en esta subrutina, la cosa no parece producir ningun beneficio)
        IMPLICIT NONE
        INTEGER NS1,NS2,NC1,NC2
C
        INCLUDE 'redlib.inc'
        INTEGER NPPMAX
        PARAMETER(NPPMAX=101)
C
        INTEGER NPP,NFIT,NDP
        INTEGER NPIX,NPIX1,NPIX2
        INTEGER II,JJ,KK
        REAL MEAN,SIGMA
        REAL A(NCMAX,NSMAX)
        REAL RSORT(NPPMAX*NPPMAX)
        REAL PERCENT1,PERCENT2
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLK1/MEAN,SIGMA
        COMMON/BLK5/NFIT,NPP,NDP
        COMMON/BLK6/PERCENT1,PERCENT2
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C ordenamos los pixels
        KK=0
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            KK=KK+1
            RSORT(KK)=A(JJ,II)
          END DO
        END DO
        NPIX=NPP*NPP
        CALL ORDENA1F(NPIX,RSORT)
C calculamos nuevo valor de MEAN y SIGMA
        NPIX1=INT(PERCENT1*REAL(NPP*NPP)/100.+0.5)
        NPIX2=INT(PERCENT2*REAL(NPP*NPP)/100.+0.5)
        NPIX=NPIX2-NPIX1+1
C
        MEAN=0.
        DO KK=NPIX1,NPIX2
          MEAN=MEAN+RSORT(KK)
        END DO
        MEAN=MEAN/REAL(NPIX)
C
        SIGMA=0.
        IF(NPIX.GT.1)THEN
          DO KK=NPIX1,NPIX2
            SIGMA=SIGMA+(MEAN-RSORT(KK))*(MEAN-RSORT(KK))
          END DO
          SIGMA=SQRT(SIGMA/REAL(NPIX-1))
        ELSE
          WRITE(*,101)'WARNING: No. of pixels employed to calculate '//
     +     'refined values'
          WRITE(*,101)'of MEAN and SIGMA = 1. Change options (#).'
        END IF
C
101     FORMAT(A)
        END
C
C******************************************************************************
C
        SUBROUTINE DISPLAYMENU
        IMPLICIT NONE
        WRITE(*,*)
        WRITE(*,101)'Use mouse to select pixels...'
        WRITE(*,*)
        WRITE(*,101)'? - This Help'
        WRITE(*,101)'x - X interpolation'
        WRITE(*,101)'y - Y interpolation'
        WRITE(*,101)'n - Number of pixels to be employed at both sides'
        WRITE(*,101)'u - Undo (recover original data)'
        WRITE(*,101)'s - Set background and foreground'
        WRITE(*,101)'# - Special OPTIONS'
        WRITE(*,101)'r - Remove all the crosses'
        WRITE(*,101)'t - Test: show current pixel value, MEAN and SIGMA'
        WRITE(*,101)'p - Plot spectrum or spatial cross section'
        WRITE(*,101)'l - Have a LOOK'
        WRITE(*,101)'c - Continue'
        WRITE(*,101)'q - Exit (RETURN) to main menu'
        WRITE(*,*)
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Dibuja una cruz sobre el pixel JJ,II
        SUBROUTINE MARCA(JJ,II)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INTEGER JJ,II
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XC,YC
        LOGICAL LCOLOR(MAX_ID_RED)
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        XC=REAL(JJ)
        YC=REAL(II)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(LCOLOR(ITERM)) CALL PGSCI(4)
          CALL PGSLW(7)
          CALL PGMOVE(XC-0.45,YC+0.45)
          CALL PGDRAW(XC+0.45,YC-0.45)
          CALL PGMOVE(XC-0.45,YC-0.45)
          CALL PGDRAW(XC+0.45,YC+0.45)
          IF(LCOLOR(ITERM)) CALL PGSCI(5)
          CALL PGSLW(3)
          CALL PGMOVE(XC-0.45,YC+0.45)
          CALL PGDRAW(XC+0.45,YC-0.45)
          CALL PGMOVE(XC-0.45,YC-0.45)
          CALL PGDRAW(XC+0.45,YC+0.45)
          IF(LCOLOR(ITERM)) CALL PGSCI(1)
          CALL PGSLW(1)
        END DO
        END
C
C******************************************************************************
C
C Detecta si alrededor del pixel J,I hay algun otro pixel que ya ha sido
C determinado como C.R. (es decir, MASK(,)=.TRUE.)
        LOGICAL FUNCTION NEIGHBOR(J,I)
        IMPLICIT NONE
        INTEGER J,I
C
        INTEGER NPPMAX
        PARAMETER(NPPMAX=101)
        INTEGER II,JJ
        INTEGER NC1,NC2,NS1,NS2
        LOGICAL MASK(NPPMAX,NPPMAX)
        COMMON/BLK4/MASK,NC1,NC2,NS1,NS2
C
        NEIGHBOR=.FALSE.
        DO II=I-1,I+1
          DO JJ=J-1,J+1
            IF((JJ.GE.NC1).AND.(JJ.LE.NC2).AND.(II.GE.NS1).AND.
     +      (II.LE.NS2))THEN
              IF((II.EQ.I).AND.(JJ.EQ.J))THEN
              ELSE
                IF(MASK(JJ-NC1+1,II-NS1+1))THEN
                  NEIGHBOR=.TRUE.
                  RETURN
                END IF
              END IF
            END IF
          END DO
        END DO
        END
C
C******************************************************************************
C
C Presenhalizamos los pixels que sobrepasen MEAN en SIGTHRESHOLD veces SIGMA.
C Hacemos una busqueda que comienza en el pixel central, y luego, a modo
C de capas de cebolla, continua la busqueda obligando a que los pixeles
C candidatos a C.R. formen un conjunto conexo.
C NPIX es el numero de pixeles del cuadrado NPP x NPP que han sido 
C testeados (sin contar el pixel original del C.R.), mientras que KK es 
C la longitud que se va recorriendo en cada perimetro cuadrado concentrico
C sin cambiar de direccion. El punto de arranque se encuentra siempre en
C la diagonal que se encuentra abajo y a la izquierda del punto de arranque.
C MASKSKY () indica si en el canal JJ el programa ha encontrado una posible
C linea de cielo, en cuyo caso no se presenhaliza ningun pixel de dicha
C columna.
C Si LPCROSS=.TRUE. dibujamos una cruz como presenhalizacion.
C Hemos anhadido la posibilidad de eliminar la forma promedio en X e Y para
C facilitar la busqueda de los rayos cosmicos. Para ello, la subrutina, cada
C que encuentra un nuevo pixel candidato a rayo cosmico, recalcula la forma
C promedio en X e/o Y e itera.
C Para calcular la forma espectral promedio (si se utiliza), eliminamos en los
C calculos el pixel con la mayor y la menor senhal.
        SUBROUTINE WALKER(J,I,LPCROSS)
        IMPLICIT NONE
        INTEGER I,J
        LOGICAL LPCROSS
C
        INCLUDE 'redlib.inc'
        INTEGER NPPMAX
        PARAMETER(NPPMAX=101)
C
        INTEGER II,JJ,KK,LL,KK0
        INTEGER NS1,NS2,NC1,NC2
        INTEGER NPP,NFIT,NDP
        INTEGER NPIX
        INTEGER CAUTOMODE
        REAL MEAN,SIGMA
        REAL A(NCMAX,NSMAX)
        REAL DEV
        REAL SIGTHRESHOLD
        REAL RQSKY
        REAL MEANX(NCMAX),MEANY(NSMAX)
        REAL RSORT(NPPMAX)
        CHARACTER*1 CCROS
        LOGICAL MASK(NPPMAX,NPPMAX),MASKSKY(NPPMAX)
        LOGICAL NEIGHBOR
        LOGICAL CRFOUND
        LOGICAL CLEANAUTO
        LOGICAL FIRSTCR
        LOGICAL INSIDE1
        LOGICAL LREMOVEX,LREMOVEY
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLK1/MEAN,SIGMA
        COMMON/BLK5/NFIT,NPP,NDP
        COMMON/BLK4/MASK,NC1,NC2,NS1,NS2
        COMMON/BLK7A/SIGTHRESHOLD
        COMMON/BLK7B/RQSKY
        COMMON/BLK8/MASKSKY
        COMMON/BLK9/CLEANAUTO,CAUTOMODE,CCROS
        COMMON/BLK18/LREMOVEX,LREMOVEY
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C FIRSTCR queda establecido a .FALSE. hasta que encontremos el primer pixel
C candidato a C.R. De esta forma la subrutina puede detectar rayos cosmicos
C que no esten centrados en la imagen de NPP x NPP pixels (lo que puede
C ocurrir si el usuario esta haciendo un zoom desde la imagen completa, y
C el cursor no estaba situado perfectamente sobre el rayo cosmico).
C inicializamos: ningun pixel es rayos cosmico
        DO II=1,NPP
          DO JJ=1,NPP
            MASK(JJ,II)=.FALSE.
          END DO
        END DO
C Pixel central
        FIRSTCR=.FALSE.
        DEV=(A(J,I)-MEAN)/SIGMA
        IF(DEV.GT.SIGTHRESHOLD)THEN
          MASK(J-NC1+1,I-NS1+1)=.TRUE.
          FIRSTCR=.TRUE.
          IF(LPCROSS) CALL MARCA(J,I)
        END IF
C forma promedio en X eliminando en el calculo los pixels que han sido
C marcados como rayo cosmico (si todo los pixels han sido marcados, se
C calcula el promedio con todos ellos) asi como los pixels en el mismo scan
C que el pixel de arranque
10      DO JJ=NC1,NC2
          MEANX(JJ)=0.
        END DO
        IF(LREMOVEX)THEN
          DO JJ=NC1,NC2
            KK=0
            DO II=NS1,NS2
              IF((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(II.NE.I))THEN
                KK=KK+1
                RSORT(KK)=A(JJ,II)
              END IF
            END DO
            IF(KK.EQ.0)THEN
              DO II=NS1,NS2
                MEANX(JJ)=MEANX(JJ)+A(JJ,II)
              END DO
              MEANX(JJ)=MEANX(JJ)/REAL(NS2-NS1+1)
            ELSEIF(KK.EQ.1)THEN
              MEANX(JJ)=RSORT(1)
            ELSEIF(KK.EQ.2)THEN
              MEANX(JJ)=(RSORT(1)+RSORT(2))/2.
            ELSE
              CALL ORDENA1F(KK,RSORT)
              DO KK0=2,KK-1
                MEANX(JJ)=MEANX(JJ)+RSORT(KK0)
              END DO
              MEANX(JJ)=MEANX(JJ)/REAL(KK-2)
            END IF
          END DO
        END IF
C forma promedio en Y eliminando en el calculo los pixels que han sido
C marcados como rayo cosmico (si todo los pixels han sido marcados, se
C calcula el promedio con todos ellos) asi como los pixels en el mismo canal
C que el pixel del arranque
        DO II=NS1,NS2
          MEANY(II)=0.
        END DO
        IF(LREMOVEY)THEN
          DO II=NS1,NS2
            KK=0
            DO JJ=NC1,NC2
              IF((.NOT.MASK(JJ-NC1+1,II-NS1+1)).AND.(JJ.NE.J))THEN
                KK=KK+1
                RSORT(KK)=A(JJ,II)-MEANX(JJ)
              END IF
            END DO
            IF(KK.EQ.0)THEN
              DO JJ=NC1,NC2
                MEANY(II)=MEANY(II)+A(JJ,II)-MEANX(JJ)
              END DO
              MEANY(II)=MEANY(II)/REAL(NC2-NC1+1)
            ELSEIF(KK.EQ.1)THEN
              MEANY(II)=RSORT(1)
            ELSEIF(KK.EQ.2)THEN
              MEANY(II)=(RSORT(1)+RSORT(2))/2.
            ELSE
              CALL ORDENA1F(KK,RSORT)
              DO KK0=2,KK-1
                MEANY(II)=MEANY(II)+RSORT(KK0)
              END DO
              MEANY(II)=MEANY(II)/REAL(KK-2)
            END IF
          END DO
        END IF
C Movimiento en coronas cuadradas de tamanho creciente. Cada vez que se 
C detecte un nuevo rayo cosmico que no hubiera sido detectado anteriormente,
C y si LREMOVEX=.TRUE. or LREMOVEY=.TRUE., se vuelve a recalcular la(s)
C forma(s) promedio correspondiente(s).
        NPIX=0
        KK=0
        CRFOUND=.TRUE.
        DO WHILE(((NPIX.LT.(NPP*NPP-1)).AND.(CRFOUND)).OR.
     +  (.NOT.FIRSTCR))
          CRFOUND=.FALSE.
          KK=KK+2
C posicion de arranque (sobre la diagonal)
          JJ=J-KK/2
          II=I-KK/2
C nos movemos hacia arriba (mientras no tengamos linea de cielo)
          DO LL=1,KK
            II=II+1
            IF(INSIDE1(JJ,NC1,NC2,II,NS1,NS2))THEN
              NPIX=NPIX+1
              IF(.NOT.MASKSKY(JJ-NC1+1))THEN
                IF(LREMOVEX.OR.LREMOVEY)THEN
                  DEV=(A(JJ,II)-MEANX(JJ)-MEANY(II))/SIGMA
                ELSE
                  DEV=(A(JJ,II)-MEAN)/SIGMA
                END IF
                IF(DEV.GT.SIGTHRESHOLD)THEN
                  IF((.NOT.FIRSTCR).OR.NEIGHBOR(JJ,II))THEN
                    IF(.NOT.MASK(JJ-NC1+1,II-NS1+1))THEN
                      MASK(JJ-NC1+1,II-NS1+1)=.TRUE.
                      IF(LPCROSS) CALL MARCA(JJ,II)
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                      IF(LREMOVEX.OR.LREMOVEY)GOTO 10
                    ELSE
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END DO
C nos movemos hacia la derecha (testeando en cada canal si hay linea de cielo)
          DO LL=1,KK
            JJ=JJ+1
            IF(INSIDE1(JJ,NC1,NC2,II,NS1,NS2))THEN
              NPIX=NPIX+1
              IF(.NOT.MASKSKY(JJ-NC1+1))THEN
                IF(LREMOVEX.OR.LREMOVEY)THEN
                  DEV=(A(JJ,II)-MEANX(JJ)-MEANY(II))/SIGMA
                ELSE
                  DEV=(A(JJ,II)-MEAN)/SIGMA
                END IF
                IF(DEV.GT.SIGTHRESHOLD)THEN
                  IF((.NOT.FIRSTCR).OR.NEIGHBOR(JJ,II))THEN
                    IF(.NOT.MASK(JJ-NC1+1,II-NS1+1))THEN
                      MASK(JJ-NC1+1,II-NS1+1)=.TRUE.
                      IF(LPCROSS) CALL MARCA(JJ,II)
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                      IF(LREMOVEX.OR.LREMOVEY)GOTO 10
                    ELSE
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END DO
C nos movemos hacia abajo (mientras no tengamos linea de cielo)
          DO LL=1,KK
            II=II-1
            IF(INSIDE1(JJ,NC1,NC2,II,NS1,NS2))THEN
              NPIX=NPIX+1
              IF(.NOT.MASKSKY(JJ-NC1+1))THEN
                IF(LREMOVEX.OR.LREMOVEY)THEN
                  DEV=(A(JJ,II)-MEANX(JJ)-MEANY(II))/SIGMA
                ELSE
                  DEV=(A(JJ,II)-MEAN)/SIGMA
                END IF
                IF(DEV.GT.SIGTHRESHOLD)THEN
                  IF((.NOT.FIRSTCR).OR.NEIGHBOR(JJ,II))THEN
                    IF(.NOT.MASK(JJ-NC1+1,II-NS1+1))THEN
                      MASK(JJ-NC1+1,II-NS1+1)=.TRUE.
                      IF(LPCROSS) CALL MARCA(JJ,II)
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                      IF(LREMOVEX.OR.LREMOVEY)GOTO 10
                    ELSE
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END DO
C nos movemos hacia la izquierda, hasta la posicion de arranque (testeando 
C en cada canal si hay linea de cielo)
          DO LL=1,KK
            JJ=JJ-1
            IF(INSIDE1(JJ,NC1,NC2,II,NS1,NS2))THEN
              NPIX=NPIX+1
              IF(.NOT.MASKSKY(JJ-NC1+1))THEN
                IF(LREMOVEX.OR.LREMOVEY)THEN
                  DEV=(A(JJ,II)-MEANX(JJ)-MEANY(II))/SIGMA
                ELSE
                  DEV=(A(JJ,II)-MEAN)/SIGMA
                END IF
                IF(DEV.GT.SIGTHRESHOLD)THEN
                  IF((.NOT.FIRSTCR).OR.NEIGHBOR(JJ,II))THEN
                    IF(.NOT.MASK(JJ-NC1+1,II-NS1+1))THEN
                      MASK(JJ-NC1+1,II-NS1+1)=.TRUE.
                      IF(LPCROSS) CALL MARCA(JJ,II)
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                      IF(LREMOVEX.OR.LREMOVEY)GOTO 10
                    ELSE
                      CRFOUND=.TRUE.
                      FIRSTCR=.TRUE.
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END DO
C
          IF((NPIX.EQ.(NPP*NPP-1)).AND.(.NOT.FIRSTCR))THEN
            FIRSTCR=.TRUE.
            WRITE(*,101)'No C.R. found!'
          END IF
        END DO
101     FORMAT(A)
        END
C
C******************************************************************************
C
C Determina si J1.LE.JJ1.LE.JJ2.LE.J2
        LOGICAL FUNCTION INSIDE2(JJ1,JJ2,J1,J2)
        IMPLICIT NONE
        INTEGER J1,J2
        INTEGER JJ1,JJ2
C
        IF((J1.LE.JJ1).AND.(JJ1.LE.JJ2).AND.(JJ2.LE.J2))THEN
          INSIDE2=.TRUE.
        ELSE
          INSIDE2=.FALSE.
        END IF
        END
C
C******************************************************************************
C
        SUBROUTINE CROPTIONS
C Definimos los parametros de busqueda de C.R.
C
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        INTEGER READILIM
        REAL READF
C
        INTEGER NAUXMAX                   !numero maximo de imagenes auxiliares
        PARAMETER (NAUXMAX=3)         !si se cambia hacerlo en todos los sitios
C
        INTEGER NWIN
        INTEGER NC,NS
        INTEGER NTYPEPLOT
        INTEGER NAUXFRAME
        INTEGER I,J
        REAL SIGCR,FACTCR,MINVCR,SIGCRAUX
        REAL SIGTHRESHOLD
        REAL AUXFRAME(NCMAX,NSMAX,NAUXMAX)
        REAL FACTOR
        REAL FIXED_FG,FIXED_BG
        CHARACTER*1 COPC,CAUX,CRXY,CCRXY,CFIX_FGBG,COTHERCR
        CHARACTER*20 CDUMMY
        CHARACTER*75 AUXFILE(NAUXMAX)
        LOGICAL LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        LOGICAL LREMOVEX,LREMOVEY
C
C Variable que preservan el valor original de la imagen que estamos limpiando
        REAL TIMEXPOS0,AIRMASS0
        CHARACTER*255 OBJECT0,FITSFILE0,COMMENT0
C
        COMMON/BLK2/NS,NC,NWIN
        COMMON/BLK3A/LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        COMMON/BLK3B/SIGCR,FACTCR,MINVCR,SIGCRAUX
        COMMON/BLK7A/SIGTHRESHOLD
        COMMON/BLK15/NTYPEPLOT
        COMMON/BLK17A/NAUXFRAME
        COMMON/BLK17B/AUXFRAME
        COMMON/BLK17C/AUXFILE
        COMMON/BLK17D/STWV,DISP
        COMMON/BLK18/LREMOVEX,LREMOVEY
        COMMON/BLK20A/CFIX_FGBG
        COMMON/BLK20B/FIXED_FG,FIXED_BG
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
C
204     WRITE(*,100)'Are you employing SIGCR  (y/n) '
        COPC(1:1)=READC('y','yn')
        IF(COPC.EQ.'y')THEN
          LSIGCR=.TRUE.
          WRITE(*,100)'SIGCR '
          WRITE(CDUMMY,*) SIGCR
          SIGCR=READF(CDUMMY)
          IF(SIGCR.LE.0.)THEN
            WRITE(*,101)'ERROR: Invalid entry for SIGCR. Try again.'
            GOTO 204
          END IF
        END IF
205     WRITE(*,100)'Are you employing FACTCR (y/n) '
        COPC(1:1)=READC('n','yn')
        IF(COPC.EQ.'y')THEN
          LFACTCR=.TRUE.
          WRITE(*,100)'FACTCR '
          WRITE(CDUMMY,*) FACTCR
          FACTCR=READF(CDUMMY)
          IF(FACTCR.LE.0.)THEN
            WRITE(*,101)'ERROR: Invalid entry for FACTCR. Try again.'
            GOTO 205
          END IF
        END IF
206     WRITE(*,100)'Are you employing MINVCR (y/n) '
        COPC(1:1)=READC('n','yn')
        IF(COPC.EQ.'y')THEN
          LMINVCR=.TRUE.
          WRITE(*,100)'MINVCR '
          WRITE(CDUMMY,*) MINVCR
          MINVCR=READF(CDUMMY)
          IF(MINVCR.LE.0.)THEN
            WRITE(*,101)'ERROR: Invalid entry for MINVCR. Try again.'
            GOTO 206
          END IF
        END IF
        IF((.NOT.LSIGCR).AND.(.NOT.LFACTCR).AND.(.NOT.LMINVCR))THEN
          WRITE(*,101)'ERROR: no method to remove C.R. has been given.'
          STOP
        END IF
C
        WRITE(*,100)'Remove mean x/y direction before looking '//
     +   'for the C.R. (y/n) '
        IF(LREMOVEX.OR.LREMOVEY)THEN
          CRXY(1:1)=READC('y','yn')
        ELSE
          CRXY(1:1)=READC('n','yn')
        END IF
        IF(CRXY.EQ.'y')THEN
          WRITE(*,100)'[x],[y] direction or [b]oth (x/y/b) '
          IF(LREMOVEX.AND.LREMOVEY)THEN
            CCRXY(1:1)=READC('b','xyb')
          ELSEIF(LREMOVEX)THEN
            CCRXY(1:1)=READC('x','xyb')
          ELSEIF(LREMOVEY)THEN
            CCRXY(1:1)=READC('y','xyb')
          ELSE
            CCRXY(1:1)=READC('b','xyb')
          END IF
          IF(CCRXY.EQ.'y')THEN
            LREMOVEX=.FALSE.
            LREMOVEY=.TRUE.
          ELSEIF(CCRXY.EQ.'x')THEN
            LREMOVEX=.TRUE.
            LREMOVEY=.FALSE.
          ELSE
            LREMOVEX=.TRUE.
            LREMOVEY=.TRUE.
          END IF
        ELSE
          LREMOVEX=.FALSE.
          LREMOVEY=.FALSE.
        END IF
C
        WRITE(*,*)
207     WRITE(*,101)'Edge size (pixels) of the Searching square'
        WRITE(*,100)'(Note: this number must be odd and .ge. 5) '
        WRITE(CDUMMY,*) NWIN
        NWIN=READI(CDUMMY)
        IF(MOD(NWIN,2).EQ.0)THEN
          NWIN=NWIN+1
          WRITE(*,100)'WARNING: Effective Edge Size changed to '
          WRITE(*,*)NWIN
        END IF
        IF(NWIN.LT.5)THEN
          WRITE(*,101)'ERROR: Invalid number (< 5). Try again.'
          GOTO 207
        END IF
        IF((NWIN.GT.NS).OR.(NWIN.GT.NC))THEN
          WRITE(*,101)'ERROR: Edge size greater than image size. '//
     +     'Try again.'
          GOTO 207
        END IF
C
208     WRITE(*,100)'SIGTHRESHOLD '
        WRITE(CDUMMY,*) SIGTHRESHOLD
        SIGTHRESHOLD=READF(CDUMMY)
        IF(SIGTHRESHOLD.LE.0.)THEN
          WRITE(*,101)'ERROR: Invalid entry for SIGTHRESHOLD. '//
     +     'Try again.'
          GOTO 208
        END IF
C
        WRITE(*,*)
        WRITE(*,101)'(1) 3D-plot BARS (hidden lines)'
        WRITE(*,101)'(2) 3D-plot BARS (transparent)'
        WRITE(*,101)'(3) image'
        WRITE(CDUMMY,*)NTYPEPLOT
        WRITE(*,100)'Option '
        NTYPEPLOT=READILIM(CDUMMY,1,3)
C
        WRITE(*,*)
        WRITE(*,100)'Are you using fixed FG and BG (y/n) '
        CFIX_FGBG(1:1)=READC(CFIX_FGBG,'yn')
        IF(CFIX_FGBG.EQ.'y')THEN
          WRITE(*,100)'Background'
          FIXED_BG=READF('@')
          WRITE(*,100)'Foreground'
          FIXED_FG=READF('@')
        END IF
C
209     WRITE(*,*)
        WRITE(*,100)'Load auxiliary frame (y/n) '
        CAUX(1:1)=READC('n','yn')
        IF(CAUX.EQ.'y')THEN
          NAUXFRAME=NAUXFRAME+1
          IF(NAUXFRAME.GT.NAUXMAX) NAUXFRAME=NAUXMAX
          WRITE(*,110)'* Next will be auxiliary frame #',NAUXFRAME
C
          OBJECT0=OBJECT
          FITSFILE0=FITSFILE
          COMMENT0=COMMENT
          TIMEXPOS0=TIMEXPOS
          AIRMASS0=AIRMASS
C
          WRITE(*,100)'Auxiliary frame file name'
          AUXFILE(NAUXFRAME)=
     +     INFILEX(20,'@',NS,NC,STWV,DISP,21,.FALSE.) !...................match
          DO I=1,NS
            READ(20) (AUXFRAME(J,I,NAUXFRAME),J=1,NC)
          END DO
C
          CLOSE(20)
          WRITE(*,100)'Factor (Auxframe=Factor*Auxframe) '
          FACTOR=READF('1.0')
          DO I=1,NS
            DO J=1,NC
              AUXFRAME(J,I,NAUXFRAME)=FACTOR*AUXFRAME(J,I,NAUXFRAME)
            END DO
          END DO
C
          OBJECT=OBJECT0
          FITSFILE=FITSFILE0
          COMMENT=COMMENT0
          TIMEXPOS=TIMEXPOS0
          AIRMASS=AIRMASS0
C
          GOTO 209
        END IF
C
        IF(NAUXFRAME.GT.1)THEN
          WRITE(*,100)'Are you activating LOTHERCR search (y/n) '
          IF(LOTHERCR)THEN
            COTHERCR='y'
          ELSE
            COTHERCR='n'
          END IF
          COTHERCR(1:1)=READC(COTHERCR,'yn')
          LOTHERCR=(COTHERCR.EQ.'y')
          IF(LOTHERCR)THEN
            WRITE(*,100)'SIGCRAUX '
210         WRITE(CDUMMY,*) SIGCRAUX
            SIGCRAUX=READF(CDUMMY)
            IF(SIGCRAUX.LE.0.)THEN
              WRITE(*,101)'ERROR: Invalid entry for SIGCRAUX. '//
     +         'Try again.'
              GOTO 210
            END IF
          END IF
        END IF
C
        RETURN
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C Permite dibujar la imagen completa o hacer zoom sobre ella
C Los parametros que se pasan a traves del COMMON son
C A(j,i) - la matriz imagen con i:scans, j:canales
C NSCAN: el numero de scans de la imagen
C NCHAN: el numero de canales de la imagen
        SUBROUTINE LOOK(FILENAME,NC1,NC2,NS1,NS2)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        INTEGER READI
        REAL READF
C
        CHARACTER*(*) FILENAME
        INTEGER NBOTONES
        INTEGER NPPMAX
        INTEGER NAUXMAX                   !numero maximo de imagenes auxiliares
        PARAMETER(NBOTONES=18)
        PARAMETER(NPPMAX=101)
        PARAMETER (NAUXMAX=3)         !si se cambia hacerlo en todos los sitios
C
        INTEGER L
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NS,NC,NWIN
        INTEGER II,JJ
        INTEGER IIMIN,IIMAX,JJMIN,JJMAX
        INTEGER IXC1,IXC2,IYC1,IYC2
        INTEGER NTPT
        INTEGER NB,NBLOCAL
        INTEGER NFIT,NPP,NDP
        INTEGER BMODE3(NBOTONES)
        INTEGER NTOTCR
        INTEGER IWIDTHPAN,JWIDTHPAN,NWPAN,NWMAX
        INTEGER I1PAN,I2PAN,J1PAN,J2PAN
        INTEGER JUSTMODE
        INTEGER BMODE2(NBOTONES)
        INTEGER NAUXFRAME,NAUX,NCTOT0,NCTOT(NAUXMAX)
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX)
        REAL AA(NCMAX,NSMAX)    !imagen backup original para [Aux.Frame] button
        REAL BG,FG,BGOLD,FGOLD
        REAL TR(6),XC,YC
        REAL MEANVAL,MAXVAL,MINVAL,SIGMAVAL
        REAL MEAN,SIGMA
        REAL SIGCR,FACTCR,MINVCR,SIGCRAUX
        REAL FACTOR
        REAL X1DUM,X2DUM,Y1DUM,Y2DUM
        REAL AUXFRAME(NCMAX,NSMAX,NAUXMAX)
        REAL X1VPORT,X2VPORT,Y1VPORT,Y2VPORT
        CHARACTER*1 CH,CBFP,CFACTORSIMUL
        CHARACTER*20 LABEL2(NBOTONES),CDUMMY
        CHARACTER*20 LABEL3(NBOTONES)
        CHARACTER*75 AUXFILE(NAUXMAX)
        CHARACTER*80 GLABEL
        LOGICAL LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        LOGICAL INSIDE1
        LOGICAL LNEXT
        LOGICAL LBEXIST
        LOGICAL LPANORAMA
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLK1/MEAN,SIGMA  
        COMMON/BLK2/NS,NC,NWIN
        COMMON/BLK3A/LSIGCR,LFACTCR,LMINVCR,LOTHERCR
        COMMON/BLK3B/SIGCR,FACTCR,MINVCR,SIGCRAUX
        COMMON/BLK5/NFIT,NPP,NDP
        COMMON/BLK10/NTOTCR
        COMMON/BLK11/LABEL2,BMODE2
        COMMON/BLK16/FACTOR
        COMMON/BLK17A/NAUXFRAME
        COMMON/BLK17B/AUXFRAME
        COMMON/BLK17C/AUXFILE
        COMMON/BLK17D/STWV,DISP
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        NB=0
        CFACTORSIMUL='n'
C
        DATA (LABEL3(II),II=1,NBOTONES)/
     +   '[z]oom (m)','zoom [k]','[w]hole','[s]et BG/FG','[p]anorama',
     +   'e[x]it','[r]egion','[j](jump)','[p](prev.)','[n](next)',
     +   '(s[t]op)','Min[,/]Max','[1] Min+','[2] Min-','[3] Max+',
     +   '[4] Max-',' ','Aux.[f]rame'/
        DATA (BMODE3(II),II=1,NBOTONES)/ 0, 0, 0, 0, 0, 0,
     +                                   0, 3, 3, 3, 3, 0,
     +                                   0, 0, 0, 0, 0, 3/
C
        WRITE(CDUMMY,*) FACTOR
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        LABEL3(17)=CDUMMY(1:L)
C
        IF(NAUXFRAME.GT.0) BMODE3(18)=0
C
        TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
        LPANORAMA=.FALSE.
        JUSTMODE=0
C
C calcula el Background y el Foreground a partir del valor medio y la
C desviacion estandard de una region alrededor del pixel central de la
C imagen
        MINVAL=A(1,1)
        MAXVAL=MINVAL
        DO II=1,NSCAN
          DO JJ=1,NCHAN
            IF(A(JJ,II).LT.MINVAL) MINVAL=A(JJ,II)
            IF(A(JJ,II).GT.MAXVAL) MAXVAL=A(JJ,II)
          END DO
        END DO
        CALL STATISTICS(NCHAN/2,NSCAN/2)
        IF(SIGMA.GT.0.0)THEN
          IF(LSIGCR)THEN
            BG=MEAN-SIGCR*SIGMA
            FG=MEAN+SIGCR*SIGMA
          ELSEIF(LFACTCR)THEN
            BG=(MEAN-1)*FACTCR
            FG=(MEAN+1)*FACTCR
          ELSEIF(LMINVCR)THEN
            BG=MEAN+MINVCR
            FG=MEAN-MINVCR
          ELSE
            BG=MEAN-SIGCRAUX*SIGMA
            FG=MEAN+SIGCRAUX*SIGMA
          END IF
          IF(BG.LT.MINVAL) BG=MINVAL
          IF(FG.GT.MAXVAL) FG=MAXVAL
        ELSE
          BG=MINVAL
          FG=MAXVAL
        END IF
C
4       NC1=1
        NC2=NCHAN
        NS1=1
        NS2=NSCAN
C recompute bf and fg using zscale
        CALL ZSCALE(A,NC1,NC2,NS1,NS2,BG,FG)
C
5       CALL PGBBUF
        CALL PGSAVE
        DO II=1,NBOTONES
          IF(BMODE3(II).GE.0) CALL BUTTON(II,LABEL3(II),0)
          IF(BMODE3(II).GE.1) CALL BUTTON(II,LABEL3(II),BMODE3(II))
        END DO
        CALL PGEBUF
        CALL PGUNSA
C
        WRITE(*,100)'Background: '
        WRITE(*,*)BG
        WRITE(*,100)'Foreground: '
        WRITE(*,*)FG
C
16      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
ccc         CALL RPGERAS
            CALL RPGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +       REAL(NS2)+.6,JUSTMODE,0)
          ELSE
            CALL PGENV(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
     +       REAL(NS2)+.6,JUSTMODE,-2)
ccc            CALL PGWINDOW(REAL(NC1)-.6,REAL(NC2)+.6,REAL(NS1)-.6,
ccc     +       REAL(NS2)+.6)
            CALL PGBOX('BCNITS',0.0,0,'BCNITS',0.0,0)
          END IF
          CALL PGMTEXT('L',2.5,.5,.5,'scan')
          CALL PGMTEXT('B',3.,.5,.5,'channel')
          CALL PGIDEN_RED
          IF(.NOT.LPANORAMA)THEN
            CALL PGMTEXT('T',1.,.5,.5,FILENAME)
          ELSE
            CALL PGMTEXT('T',1.,0.,0.,FILENAME)
            WRITE(GLABEL,'(I6,A,I6)')NWPAN,'/',NWMAX
            CALL RMBLANK(GLABEL,GLABEL,L)
            CALL PGMTEXT('T',1.,1.,1.,GLABEL(1:L))
          END IF
        END DO
        MINVAL=A(NC1,NS1)
        MAXVAL=MINVAL
        IIMIN=NS1
        IIMAX=NS1
        JJMIN=NC1
        JJMAX=NC1
        MEANVAL=0.
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
            MEANVAL=MEANVAL+A(JJ,II)
          END DO
        END DO
        NTPT=(NS2-NS1+1)*(NC2-NC1+1)
        MEANVAL=MEANVAL/REAL((NS2-NS1+1)*(NC2-NC1+1))
        WRITE(*,'(A,I5,A,I5)')'> From Scan    #',NS1,' to ',NS2
        WRITE(*,'(A,I5,A,I5)')'> From Channel #',NC1,' to ',NC2
        WRITE(*,'(A,I10)')'> Total number of pixels: ',NTPT 
        WRITE(*,100)'> Maximum: '
        WRITE(CDUMMY,*)MAXVAL
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)CDUMMY(1:L)
        WRITE(*,100)'  at (channel,scan):'
        WRITE(*,*)JJMAX,IIMAX
        WRITE(*,100)'> Minimum: '
        WRITE(CDUMMY,*)MINVAL
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)CDUMMY(1:L)
        WRITE(*,100)'  at (channel,scan):'
        WRITE(*,*)JJMIN,IIMIN
        WRITE(*,100)'> Mean   : '
        WRITE(*,*)MEANVAL
        SIGMAVAL=0. 
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            SIGMAVAL=SIGMAVAL+(A(JJ,II)-MEANVAL)*(A(JJ,II)-MEANVAL) 
          END DO
        END DO
        SIGMAVAL=SQRT(SIGMAVAL/REAL(NTPT))
        WRITE(*,100)'> Sigma  : '
        WRITE(*,*)SIGMAVAL
        WRITE(*,*)
        IF(CBFP.NE.'y')THEN
          IF(LPANORAMA)THEN
            BG=MINVAL
            FG=MAXVAL
          END IF
        END IF
17      DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          !marcamos los puntos con el valor máximo y el mínimo
          CALL PGSCI(3)
          CALL PGPOINT(1,REAL(JJMAX),REAL(IIMAX),23)
          CALL PGSCI(5)
          CALL PGPOINT(1,REAL(JJMIN),REAL(IIMIN),23)
          CALL PGSCI(1)
        END DO
        IF(.NOT.LPANORAMA)THEN
          WRITE(CDUMMY,*)BG
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
            DO WHILE(CDUMMY(L:L).EQ.'0')
              L=L-1
            END DO
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
          END DO
          WRITE(CDUMMY,*)FG
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
            DO WHILE(CDUMMY(L:L).EQ.'0')
              L=L-1
            END DO
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
          END DO
        END IF
C------------------------------------------------------------------------------
C
18      CONTINUE
C------------------------------------------------------------------------------
20      NB=0
        CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
        CALL IFBUTTON(XC,YC,NB)
        IF(LPANORAMA)THEN
          IF(((NB.LT.8).OR.(NB.GT.11)).AND.(NB.NE.4))THEN
            IF(CH.EQ.'j')THEN
              NB=8
            ELSEIF(CH.EQ.'p')THEN
              NB=9
            ELSEIF(CH.EQ.'n')THEN
              NB=10
            ELSEIF(CH.EQ.'t')THEN
              NB=11
            END IF
          END IF
        END IF
C
        NBLOCAL=INDEX('zkwspxr    ,1234 f',CH)
        IF((NBLOCAL.EQ.0).AND.(CH.EQ.'/')) NBLOCAL=12  !zmin,zmax
        IF((NBLOCAL.NE.0).AND.(CH.NE.' '))THEN
          CALL BUTTQEX(NBLOCAL,LBEXIST)
          IF(LBEXIST) NB=NBLOCAL
        END IF
C------------------------------------------------------------------------------
        IF(NB.EQ.1)THEN
          CALL BUTTON(1,LABEL3(1),5)
          WRITE(*,101)'Press cursor at two corners of the imaginary '
     +     //'BOX to be zoomed'
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(1,LABEL3(1),0)
            GOTO 18
          END IF
          IXC1=INT(XC+0.5)
          IYC1=INT(YC+0.5)
          IF(IXC1.LT.NC1) IXC1=NC1
          IF(IYC1.LT.NS1) IYC1=NS1
          IF(IXC1.GT.NC2) IXC1=NC2
          IF(IYC1.GT.NS2) IYC1=NS2
          WRITE(*,117)'Cursor at ',IXC1,IYC1
C
          IF(LCOLOR(1)) CALL PGSCI(2)
          CALL RPGBAND(2,0,REAL(IXC1),REAL(IYC1),XC,YC,CH)
          IF(LCOLOR(1)) CALL PGSCI(1)
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            CALL BUTTON(1,LABEL3(1),0)
            GOTO 18
          END IF
          IXC2=INT(XC+0.5)
          IYC2=INT(YC+0.5)
          IF(IXC2.LT.NC1) IXC2=NC1
          IF(IYC2.LT.NS1) IYC2=NS1
          IF(IXC2.GT.NC2) IXC2=NC2
          IF(IYC2.GT.NS2) IYC2=NS2
          WRITE(*,117)'Cursor at ',IXC2,IYC2
C
          NC1=MIN0(IXC1,IXC2)
          NC2=MAX0(IXC1,IXC2)
          NS1=MIN0(IYC1,IYC2)
          NS2=MAX0(IYC1,IYC2)
          CALL BUTTON(1,LABEL3(1),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.2)THEN
          CALL BUTTON(2,LABEL3(2),5) 
          WRITE(*,101)'Enter the coordinates of two corners of the '//
     +     'imaginary box to be zoomed.'
          WRITE(*,100)'First point  (channel,scan) '
          CALL READ2I('1,1',IXC1,IYC1)
          IF(.NOT.INSIDE1(IXC1,1,NCHAN,IYC1,1,NSCAN))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(2,LABEL3(2),0) 
            GOTO 18
          END IF
          WRITE(*,100)'Second point (channel,scan) '
          WRITE(CDUMMY,'(I9,A1,I9)')NCHAN,',',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL READ2I(CDUMMY(1:L),IXC2,IYC2)
          IF(.NOT.INSIDE1(IXC2,1,NCHAN,IYC2,1,NSCAN))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(2,LABEL3(2),0) 
            GOTO 18
          END IF
          NC1=MIN0(IXC1,IXC2)
          NC2=MAX0(IXC1,IXC2)
          NS1=MIN0(IYC1,IYC2)
          NS2=MAX0(IYC1,IYC2)
          CALL BUTTON(2,LABEL3(2),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.3)THEN
          CALL BUTTON(3,LABEL3(3),5)
          NC1=1
          NC2=NCHAN
          NS1=1
          NS2=NSCAN
          CALL BUTTON(3,LABEL3(3),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.4)THEN
          CALL BUTTON(4,LABEL3(4),5)
          BGOLD=BG
          FGOLD=FG
          WRITE(*,100)'Background    : '
          WRITE(*,*)BG
          WRITE(*,100)'Foreground    : '
          WRITE(*,*)FG
          WRITE(*,100)'New background '
          WRITE(CDUMMY,*)BG
          BG=READF(CDUMMY)
          WRITE(*,100)'New foreground '
          WRITE(CDUMMY,*)FG
          FG=READF(CDUMMY)
          IF(.NOT.LPANORAMA)THEN
            WRITE(CDUMMY,*)BGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
            END DO
            WRITE(CDUMMY,*)FGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
              CALL PGSCI(1)
            END DO
          END IF
          CALL BUTTON(4,LABEL3(4),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.5)THEN
          LPANORAMA=.TRUE.
          JUSTMODE=1
          CALL BUTTON(8,LABEL3(8),0) 
          CALL BUTTON(9,LABEL3(9),0) 
          CALL BUTTON(10,LABEL3(10),0)
          CALL BUTTON(11,LABEL3(11),0)
          DO II=1,NBOTONES
            IF(II.LT.8)THEN
              IF(II.NE.4) CALL BUTTON(II,LABEL3(II),3)
            END IF
          END DO
209       WRITE(*,100)'Width of the panorama window (pixels) '
          IF(NSCAN.LT.NCHAN)THEN
            WRITE(CDUMMY,*) NSCAN
          ELSE
            WRITE(CDUMMY,*) NCHAN
          END IF
          IWIDTHPAN=READI(CDUMMY)
          IF((IWIDTHPAN.LT.1).OR.(IWIDTHPAN.GT.NS).OR.
     +     (IWIDTHPAN.GT.NC))THEN
            WRITE(*,101)'ERROR: number out of range. Try again.'
            GOTO 209
          END IF
          CALL BUTTQPR(X1VPORT,X2VPORT,Y1VPORT,Y2VPORT)
ccc       CALL PGSLCT(IDN(1))
          CALL PGVPORT(X1VPORT,X2VPORT,Y1VPORT,Y2VPORT)
          CALL PGQVP(1,X1DUM,X2DUM,Y1DUM,Y2DUM)
          JWIDTHPAN=NINT(REAL(IWIDTHPAN)*
     +     ABS((X2DUM-X1DUM)/(Y2DUM-Y1DUM)))
          IF(JWIDTHPAN.GT.NC) JWIDTHPAN=NC
          WRITE(*,100)'Fixed BG/FG (y/n) '
          CBFP(1:1)=READC('y','yn')
          NWPAN=1
          CALL WHICHPAN(NWPAN,IWIDTHPAN,JWIDTHPAN,J1PAN,I1PAN,J2PAN,
     +     I2PAN,NWMAX)
          IF(NWMAX.GT.999)THEN
            WRITE(*,101)'ERROR: width of panorama window too small.'
            WRITE(*,101)'Try again.'
            GOTO 209
          END IF
          WRITE(*,'(A,I3,A,I3)')'Current panorama window ',NWPAN,
     +     ' of ',NWMAX
          NC1=J1PAN
          NC2=J2PAN
          NS1=I1PAN
          NS2=I2PAN
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.6)THEN
          CALL BUTTON(6,LABEL3(6),5)
          RETURN
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.7)THEN
          CALL BUTTON(7,LABEL3(7),5)
          WRITE(*,100)'Channel, Scan to go'
          CALL READ2I('@',IXC1,IYC1)
          IF(.NOT.INSIDE1(IXC1,NC1,NC2,IYC1,NS1,NS2))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            CALL BUTTON(7,LABEL3(7),0)
            GOTO 18
          END IF
          CALL BUTTON(7,LABEL3(7),0)
          CALL STATISTICS(IXC1,IYC1)
          CALL PGBBUF
          CALL PGSAVE
          CALL PGERAS
ccc       CALL RPGERAS
ccc       CALL RPGERASB
          DO II=1,NBOTONES
            CDUMMY=LABEL2(II)
            IF(II.EQ.6)THEN
              WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL BUTTON(II,CDUMMY(1:L),0)
            ELSE
              IF(BMODE2(II).GE.0) CALL BUTTON(II,CDUMMY,0)
              IF(BMODE2(II).GT.0) CALL BUTTON(II,CDUMMY,BMODE2(II))
            END IF
          END DO
          CALL PGEBUF
          CALL PGUNSA
          CALL SUBPLOT(IXC1,IYC1,LNEXT,.TRUE.)
          CALL RPGERASB
          CALL BUTTSPR(0.10,0.95,0.10,0.70)
          GOTO 5
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.8)THEN
          CALL BUTTON(8,LABEL3(8),5) 
          WRITE(*,'(A,I3,A,I3)')'Current panorama window ',NWPAN,
     +     ' of ',NWMAX
210       WRITE(*,100)'New window number'
          NWPAN=READI('@')
          IF((NWPAN.LT.1).OR.(NWPAN.GT.NWMAX))THEN
            WRITE(*,101)'ERROR: invalid window number. Try again.'
            GOTO 210
          END IF
          CALL WHICHPAN(NWPAN,IWIDTHPAN,JWIDTHPAN,J1PAN,I1PAN,J2PAN,
     +     I2PAN,NWMAX)
          NC1=J1PAN
          NC2=J2PAN
          NS1=I1PAN
          NS2=I2PAN
          CALL BUTTON(8,LABEL3(8),0) 
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.9)THEN
          CALL BUTTON(9,LABEL3(9),5) 
          NWPAN=NWPAN-1
          IF(NWPAN.LT.1) NWPAN=NWMAX
          CALL WHICHPAN(NWPAN,IWIDTHPAN,JWIDTHPAN,J1PAN,I1PAN,J2PAN,
     +     I2PAN,NWMAX)
          WRITE(*,'(A,I3,A,I3)')'Current panorama window ',NWPAN,
     +     ' of ',NWMAX
          NC1=J1PAN
          NC2=J2PAN
          NS1=I1PAN
          NS2=I2PAN
          CALL BUTTON(9,LABEL3(9),0) 
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.10)THEN
          CALL BUTTON(10,LABEL3(10),5)
          NWPAN=NWPAN+1
          IF(NWPAN.GT.NWMAX) NWPAN=1
          CALL WHICHPAN(NWPAN,IWIDTHPAN,JWIDTHPAN,J1PAN,I1PAN,J2PAN,
     +     I2PAN,NWMAX)
          WRITE(*,'(A,I3,A,I3)')'Current panorama window ',NWPAN,
     +     ' of ',NWMAX
          NC1=J1PAN
          NC2=J2PAN
          NS1=I1PAN
          NS2=I2PAN
          CALL BUTTON(10,LABEL3(10),0)
          GOTO 16
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.11)THEN
          CALL BUTTON(11,LABEL3(11),5)
          LPANORAMA=.FALSE.
          JUSTMODE=0
          CALL BUTTON(11,LABEL3(11),0)
          CALL RPGERASB
          GOTO 4 
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.12)THEN
          CALL BUTTON(12,LABEL3(12),5)
          BGOLD=BG
          FGOLD=FG
          IF(CH.EQ.'/')THEN
            CALL ZSCALE(A,NC1,NC2,NS1,NS2,BG,FG)
          ELSE
            BG=MINVAL
            FG=MAXVAL
          END IF
          IF(.NOT.LPANORAMA)THEN
            WRITE(CDUMMY,*)BGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
            END DO
            WRITE(CDUMMY,*)FGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
              CALL PGSCI(1)
            END DO
          END IF
          CALL BUTTON(12,LABEL3(12),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.13)THEN
          CALL BUTTON(13,LABEL3(13),5)
          BGOLD=BG
          BG=BG+FACTOR
          IF(.NOT.LPANORAMA)THEN
            WRITE(CDUMMY,*)BGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
              CALL PGSCI(1)
            END DO
          END IF
C
          IF(CFACTORSIMUL.EQ.'y')THEN
            CALL BUTTON(15,LABEL3(15),5)
            FGOLD=FG
            FG=FG+FACTOR
            IF(.NOT.LPANORAMA)THEN
              WRITE(CDUMMY,*)FGOLD
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
                DO WHILE(CDUMMY(L:L).EQ.'0')
                  L=L-1
                END DO
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSCI(0)
                CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
                CALL PGSCI(1)
              END DO
            END IF
            CALL BUTTON(15,LABEL3(15),0)
          END IF
C
          CALL BUTTON(13,LABEL3(13),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.14)THEN
          CALL BUTTON(14,LABEL3(14),5)
          BGOLD=BG
          BG=BG-FACTOR
          IF(.NOT.LPANORAMA)THEN
            WRITE(CDUMMY,*)BGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
              CALL PGSCI(1)
            END DO
          END IF
C
          IF(CFACTORSIMUL.EQ.'y')THEN
            CALL BUTTON(16,LABEL3(16),5)
            FGOLD=FG
            FG=FG-FACTOR
            IF(.NOT.LPANORAMA)THEN
              WRITE(CDUMMY,*)FGOLD
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
                DO WHILE(CDUMMY(L:L).EQ.'0')
                  L=L-1
                END DO
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSCI(0)
                CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
                CALL PGSCI(1)
              END DO
            END IF
            CALL BUTTON(16,LABEL3(16),0)
          END IF
C
          CALL BUTTON(14,LABEL3(14),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.15)THEN
          CALL BUTTON(15,LABEL3(15),5)
          FGOLD=FG
          FG=FG+FACTOR
          IF(.NOT.LPANORAMA)THEN
            WRITE(CDUMMY,*)FGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
              CALL PGSCI(1)
            END DO
          END IF
C
          IF(CFACTORSIMUL.EQ.'y')THEN
            CALL BUTTON(13,LABEL3(13),5)
            BGOLD=BG
            BG=BG+FACTOR
            IF(.NOT.LPANORAMA)THEN
              WRITE(CDUMMY,*)BGOLD
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
                DO WHILE(CDUMMY(L:L).EQ.'0')
                  L=L-1
                END DO
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSCI(0)
                CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
                CALL PGSCI(1)
              END DO
            END IF
            CALL BUTTON(13,LABEL3(13),0)
          END IF
C
          CALL BUTTON(15,LABEL3(15),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.16)THEN
          CALL BUTTON(16,LABEL3(16),5)
          FGOLD=FG
          FG=FG-FACTOR
          IF(.NOT.LPANORAMA)THEN
            WRITE(CDUMMY,*)FGOLD
            CALL RMBLANK(CDUMMY,CDUMMY,L)
            IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
              DO WHILE(CDUMMY(L:L).EQ.'0')
                L=L-1
              END DO
            END IF
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGSCI(0)
              CALL PGMTEXT('T',1.,1.,1.,'FG='//CDUMMY(1:L))
              CALL PGSCI(1)
            END DO
          END IF
C
          IF(CFACTORSIMUL.EQ.'y')THEN
            CALL BUTTON(14,LABEL3(14),5)
            BGOLD=BG
            BG=BG-FACTOR
            IF(.NOT.LPANORAMA)THEN
              WRITE(CDUMMY,*)BGOLD
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
                DO WHILE(CDUMMY(L:L).EQ.'0')
                  L=L-1
                END DO
              END IF
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                CALL PGSCI(0)
                CALL PGMTEXT('T',1.,0.,0.,'BG='//CDUMMY(1:L))
                CALL PGSCI(1)
              END DO
            END IF
            CALL BUTTON(14,LABEL3(14),0)
          END IF
C
          CALL BUTTON(16,LABEL3(16),0)
          GOTO 17
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.17)THEN
          CALL BUTTON(17,LABEL3(17),5)
          WRITE(*,100)'New factor '
          FACTOR=READF(LABEL3(17))
          WRITE(CDUMMY,*) FACTOR
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          LABEL3(17)=CDUMMY(1:L)
          WRITE(*,100)'Change Max. and Min. simultaneously (y/n) '
          CFACTORSIMUL(1:1)=READC(CFACTORSIMUL,'yn')
          CALL BUTTON(17,LABEL3(17),0)
          GOTO 20
C------------------------------------------------------------------------------
        ELSEIF(NB.EQ.18)THEN
          CALL BUTTON(18,LABEL3(18),1)
          CALL PGBBUF
          CALL PGSAVE
          DO II=1,NBOTONES
            IF(II.NE.18)THEN
              CALL BUTTON(II,LABEL3(II),3)
            END IF
          END DO
          CALL PGEBUF
          CALL PGUNSA
C guardamos imagen inicial
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              AA(JJ,II)=A(JJ,II)
            END DO
          END DO
C
          NCTOT0=0
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              NCTOT0=NCTOT0+AA(JJ,II)
            END DO
          END DO
          DO NAUX=1,NAUXFRAME
            NCTOT(NAUX)=0
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                NCTOT(NAUX)=NCTOT(NAUX)+AUXFRAME(JJ,II,NAUX)
              END DO
            END DO
          END DO
C
          WRITE(*,*)
          WRITE(*,101)'* Press any key to change frame'
          WRITE(*,101)'* Press <q>, <x> or <right mouse button> to exit'
          WRITE(*,*)
          NAUX=NAUXFRAME
212       NAUX=NAUX+1
          IF(NAUX.EQ.NAUXFRAME+1)THEN
            WRITE(*,101)'* Displayed image is Original Frame'
            WRITE(*,100)'- '
            WRITE(*,101)FILENAME(1:TRUELEN(FILENAME))
            DO II=NS1,NS2
              DO JJ=NC1,NC2
                A(JJ,II)=AA(JJ,II)
              END DO
            END DO
            NAUX=0
          ELSE
            WRITE(*,110)'* Displayed image is auxiliary frame #',NAUX
            WRITE(*,101)'- '//AUXFILE(NAUX)
            DO II=NS1,NS2
              DO JJ=NC1,NC2
ccc                A(JJ,II)=AUXFRAME(JJ,II,NAUX)*NCTOT0/NCTOT(NAUX)
                A(JJ,II)=AUXFRAME(JJ,II,NAUX)
              END DO
            END DO
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
          END DO
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL CHLOWER(CH)
          IF((CH.NE.'q').AND.(CH.NE.'x')) GOTO 212
C recuperamos imagen original y la dibujamos
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              A(JJ,II)=AA(JJ,II)
            END DO
          END DO
          IF(NAUX.NE.0)THEN
            WRITE(*,101)'* Displayed image is Original Frame'
            WRITE(*,100)'- '
            WRITE(*,101)FILENAME(1:TRUELEN(FILENAME))
             DO ITERM=NTERM,1,-1
               CALL PGSLCT(IDN(ITERM))
               CALL PGIMAG(A,NCMAX,NSMAX,NC1,NC2,NS1,NS2,FG,BG,TR)
             END DO
          END IF
C
          CALL PGBBUF
          CALL PGSAVE
          IF(LPANORAMA)THEN
            CALL BUTTON(4,LABEL3(4),0)
            CALL BUTTON(8,LABEL3(8),0) 
            CALL BUTTON(9,LABEL3(9),0) 
            CALL BUTTON(10,LABEL3(10),0)
            CALL BUTTON(11,LABEL3(11),0)
          ELSE
            DO II=1,NBOTONES
              IF(II.NE.18)THEN
                IF(BMODE3(II).GE.0)CALL BUTTON(II,LABEL3(II),0)
                IF(BMODE3(II).GE.1)
     +           CALL BUTTON(II,LABEL3(II),BMODE3(II))
              END IF
            END DO
          END IF
          CALL PGEBUF
          CALL PGUNSA
          CALL BUTTON(18,LABEL3(18),0)
          GOTO 20
C------------------------------------------------------------------------------
C look for C.R. y cursor region
        ELSE
          IF((CH.NE.'A').AND.(CH.NE.'D').AND.(CH.NE.'X')) THEN
            WRITE(*,101)'ERROR: mouse buttom has not been detected.'
            GOTO 20
          END IF
          IXC1=INT(XC+0.5)
          IYC1=INT(YC+0.5)
          IF(.NOT.INSIDE1(IXC1,NC1,NC2,IYC1,NS1,NS2))THEN
            WRITE(*,101)'ERROR: limits out of plot.'
            GOTO 18
          END IF
          WRITE(*,117)'Cursor at ',IXC1,IYC1
          CALL STATISTICS(IXC1,IYC1)
          CALL PGBBUF
          CALL PGSAVE
          CALL PGERAS
ccc       CALL RPGERAS
ccc       CALL RPGERASB
          DO II=1,NBOTONES
            CDUMMY=LABEL2(II)
            IF(II.EQ.6)THEN
              WRITE(CDUMMY,'(A,I2,A,I2)')'[n]=',NFIT,',d=',NDP-1
              CALL RMBLANK(CDUMMY,CDUMMY,L)
              CALL BUTTON(II,CDUMMY(1:L),0)
            ELSE
              IF(BMODE2(II).GE.0) CALL BUTTON(II,CDUMMY,0)
              IF(BMODE2(II).GT.0) CALL BUTTON(II,CDUMMY,BMODE2(II))
            END IF
          END DO
          CALL PGEBUF
          CALL PGUNSA
          CALL SUBPLOT(IXC1,IYC1,LNEXT,.TRUE.)
          CALL RPGERASB
          CALL BUTTSPR(0.10,0.95,0.10,0.70)
          IF(LPANORAMA)THEN
            CALL BUTTON(4,LABEL3(4),0)
            CALL BUTTON(8,LABEL3(8),0)
            CALL BUTTON(9,LABEL3(9),0)
            CALL BUTTON(10,LABEL3(10),0)
            CALL BUTTON(11,LABEL3(11),0)
            CALL BUTTON(12,LABEL3(12),0)
            DO II=1,NBOTONES
              IF((II.LT.8).AND.(II.NE.4))THEN
                CALL BUTTON(II,LABEL3(II),0)
                CALL BUTTON(II,LABEL3(II),3)
              END IF
              IF(II.GT.12)THEN
                CALL BUTTON(II,LABEL3(II),0)
                IF(II.EQ.18) CALL BUTTON(II,LABEL3(II),BMODE3(II))
              END IF
            END DO
            GOTO 16
          ELSE
            GOTO 5
          END IF
        END IF
C------------------------------------------------------------------------------
        GOTO 20
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
117     FORMAT(A,I5,2X,I5)
C
        END
C
C******************************************************************************
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
C determina las esquinas (J1PAN,I1PAN) y (J2PAN,I2PAN) de la ventana 
C de panorama numero NWPAN, asi como el numero maximo de ventanas de
C panorama admisibles NWMAX, siendo el ancho de estas ventanas igual
C a IWIDTHPAN (en pixels)
        SUBROUTINE WHICHPAN(NWPAN,IWIDTHPAN,JWIDTHPAN,J1PAN,I1PAN,
     +   J2PAN,I2PAN,NWMAX)
        IMPLICIT NONE
        INTEGER NWPAN,IWIDTHPAN,JWIDTHPAN
        INTEGER I1PAN,J1PAN,I2PAN,J2PAN
        INTEGER NWMAX
C
        INTEGER NS,NC,NWIN
        INTEGER NWINX,NWINY
        COMMON/BLK2/NS,NC,NWIN
C
        NWINX=NC/JWIDTHPAN
        IF(MOD(NC,JWIDTHPAN).NE.0) NWINX=NWINX+1
        NWINY=NS/IWIDTHPAN
        IF(MOD(NS,IWIDTHPAN).NE.0) NWINY=NWINY+1
        NWMAX=NWINX*NWINY
C
        J1PAN=NWPAN-((NWPAN-1)/NWINX)*NWINX
        J1PAN=(J1PAN-1)*JWIDTHPAN
        J2PAN=J1PAN+JWIDTHPAN
        IF(J2PAN.GT.NC)THEN
          J2PAN=NC
          J1PAN=J2PAN-JWIDTHPAN+1
        ELSE
          J1PAN=J1PAN+1
        END IF
C
        I1PAN=((NWPAN-1)/NWINX)*IWIDTHPAN
        I2PAN=I1PAN+IWIDTHPAN
        IF(I2PAN.GT.NS)THEN
          I2PAN=NS
          I1PAN=I2PAN-IWIDTHPAN+1
        ELSE
          I1PAN=I1PAN+1
        END IF
        END
C
C******************************************************************************
C Realiza un dibujo tridimensional de la region de busqueda de rayo cosmico
C IMODE=0: lineas ocultas
C IMODE=1: lineas transparentes
C ANGLE puede ser 0, 90, 180 o 270 (angulo de vision)
C MM3D=.FALSE. limites automaticos (se devuelven como output en BG3D,FG3D)
C MM3D=.TRUE. limites fijos (se toman BG3D y FG3D como input)
        SUBROUTINE 
     +   PLOT3DBARS(NCC1,NCC2,NSS1,NSS2,IMODE,ANGLE,MM3D,BG3D,FG3D)
        IMPLICIT NONE
        INTEGER NCC1,NCC2,NSS1,NSS2
        INTEGER IMODE,ANGLE
        LOGICAL MM3D                           !si .FALSE., limites automaticos
        REAL BG3D,FG3D
C
        INCLUDE 'redlib.inc'
C
        INTEGER NPPMAX
        PARAMETER (NPPMAX=101)
C
        INTEGER NC1,NC2,NS1,NS2
        INTEGER I,J
        INTEGER K,L
        INTEGER NCOLOR(10)
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL A(NCMAX,NSMAX)
        REAL AA(NPPMAX,NPPMAX)
        REAL XMIN,XMAX,YMIN,YMAX,YMIN0,YMAX0,DX
        REAL XOFFSET,YOFFSET
        REAL X(5),Y(5),YB(4)
        REAL XX(7),YY(7)
        REAL F
        REAL X1VPORT,X2VPORT,Y1VPORT,Y2VPORT
        CHARACTER*50 CDUMMY
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        DATA (NCOLOR(K),K=1,10) /14,15,5,4,6,11,7,8,2,3/
C
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF((NCC2-NCC1+1).GT.NPPMAX)THEN
          WRITE(*,101) 'FATAL ERROR in subroutine PLOT3DBARS:'
          WRITE(*,100) 'NCC1,NCC2,NPPMAX: '
          WRITE(*,*) NCC1,NCC2,NPPMAX
          STOP
        END IF
        IF((NSS2-NSS1+1).GT.NPPMAX)THEN
          WRITE(*,101) 'FATAL ERROR in subroutine PLOT3DBARS:'
          WRITE(*,100) 'NSS1,NSS2,NPPMAX: '
          WRITE(*,*) NSS1,NSS2,NPPMAX
          STOP
        END IF
        IF(ANGLE.EQ.90)THEN
          NC1=NSS1
          NC2=NSS2
          NS1=NCC1
          NS2=NCC2
          DO I=NS1,NS2
            DO J=NC1,NC2
              AA(J-NC1+1,I-NS1+1)=A(I,NSS2-(J-NSS1))
            END DO
          END DO
        ELSEIF(ANGLE.EQ.180)THEN
          NC1=NCC1
          NC2=NCC2
          NS1=NSS1
          NS2=NSS2
          DO I=NS1,NS2
            DO J=NC1,NC2
              AA(J-NC1+1,I-NS1+1)=A(NCC2-(J-NCC1),NSS2-(I-NSS1))
            END DO
          END DO
        ELSEIF(ANGLE.EQ.270)THEN
          NC1=NSS1
          NC2=NSS2
          NS1=NCC1
          NS2=NCC2
          DO I=NS1,NS2
            DO J=NC1,NC2
              AA(J-NC1+1,I-NS1+1)=A(NCC2-(I-NCC1),J)
            END DO
          END DO
        ELSE                                                           !ANGLE=0
          NC1=NCC1
          NC2=NCC2
          NS1=NSS1
          NS2=NSS2
          DO I=NS1,NS2
            DO J=NC1,NC2
              AA(J-NC1+1,I-NS1+1)=A(J,I)
            END DO
          END DO
        END IF
C------------------------------------------------------------------------------
        YMIN0=AA(1,1)
        YMAX0=YMIN0
        DO I=NS1,NS2
          DO J=NC1,NC2
            IF(AA(J-NC1+1,I-NS1+1).LT.YMIN0) YMIN0=AA(J-NC1+1,I-NS1+1)
            IF(AA(J-NC1+1,I-NS1+1).GT.YMAX0) YMAX0=AA(J-NC1+1,I-NS1+1)
          END DO
        END DO
        IF(.NOT.MM3D)THEN
          BG3D=YMIN0
          FG3D=YMAX0
          IF(BG3D.EQ.FG3D)THEN
            IF(BG3D.EQ.0.)THEN
              BG3D=-1.
              FG3D=1.
            ELSE
              BG3D=0.9*BG3D
              FG3D=1.1*FG3D
            END IF
          END IF
        END IF
C
        XOFFSET=REAL(NC2-NC1+1)/(4*REAL(NS2-NS1+1))
        YOFFSET=0.8*(FG3D-BG3D)/REAL(NS2-NS1+1)
        XMIN=REAL(NC1)-1.
        XMAX=REAL(NC2)+REAL(NS2-NS1)*XOFFSET+1.
        YMAX=FG3D+REAL(NS2-NS1-1)*YOFFSET
        YMAX=YMAX+YOFFSET
        YMIN=BG3D-YOFFSET
        DX=XMAX-XMIN
C
        X1VPORT=0.10
        X2VPORT=0.50
        Y1VPORT=0.15
        Y2VPORT=0.75
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGBBUF
          CALL PGSAVE
          CALL PGVPORT(0.05,X2VPORT,0.05,Y2VPORT)
          CALL PGWINDOW(0.,1.,0.,1.)
          CALL PGSFS(1)
          CALL PGSCI(0)
          CALL PGRECT(0.,1.,0.,1.)
          CALL PGSCI(1)
          CALL PGVPORT(X1VPORT,X2VPORT,Y1VPORT,Y2VPORT)
          CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
          CALL PGSFS(1)
C------------------------------------------------------------------------------
          DO I=NS2,NS1,-1
            DO J=NC1,NC2
              X(1)=REAL(J)-.5+(REAL(I-NS1)-.5)*XOFFSET
              X(2)=X(1)+1.
              X(3)=REAL(J)+.5+(REAL(I-NS1)+.5)*XOFFSET
              X(4)=X(3)-1.
              Y(1)=AA(J-NC1+1,I-NS1+1)+(REAL(I-NS1)-.5)*YOFFSET
              Y(2)=Y(1)
              Y(3)=AA(J-NC1+1,I-NS1+1)+(REAL(I-NS1)+.5)*YOFFSET
              Y(4)=Y(3)
              X(5)=X(1)
              Y(5)=Y(1)
              YB(1)=BG3D+(REAL(I-NS1)-.5)*YOFFSET
              YB(2)=YB(1)
              YB(3)=BG3D+(REAL(I-NS1)+.5)*YOFFSET
              YB(4)=YB(3)
              F=(AA(J-NC1+1,I-NS1+1)-YMIN0)/(YMAX0-YMIN0)
              IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR(NINT(F*9)+1))
C..............................................................................
              IF(IMODE.EQ.0)THEN
                XX(1)=X(1)
                XX(2)=X(1)
                XX(3)=X(4)
                XX(4)=X(3)
                XX(5)=X(3)
                XX(6)=X(2)
                XX(7)=XX(1)
                YY(1)=YB(1)
                YY(2)=Y(1)
                YY(3)=Y(4)
                YY(4)=Y(3)
                YY(5)=YB(3)
                YY(6)=YB(2)
                YY(7)=YY(1)
                CALL PGPOLY(7,XX,YY)
                CALL PGSCI(0)
                CALL PGLINE(5,X,Y)
                DO K=1,3
                  CALL PGMOVE(X(K),Y(K))
                  CALL PGDRAW(X(K),YB(K))
                END DO
C..............................................................................
              ELSEIF(IMODE.EQ.1)THEN
                CALL PGPOLY(5,X,Y)
                DO K=1,4
                  CALL PGMOVE(X(K),Y(K))
                  CALL PGDRAW(X(K),YB(K))
                END DO
C..............................................................................
              ELSE
              END IF
            END DO
          END DO
          CALL PGEBUF
          CALL PGUNSA
        END DO
C------------------------------------------------------------------------------
        WRITE(CDUMMY,*)ANGLE
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(1)
          CALL PGMTEXT('T',-1.1,0.,0.,CDUMMY(1:L)//'\uo')
          CALL PGVPORT(X1VPORT,X2VPORT,Y1VPORT-0.04,Y1VPORT-0.02)
          CALL PGWINDOW(0.,10.,0.,1.)
          DO K=1,10
            IF(LCOLOR(ITERM)) CALL PGSCI(NCOLOR(K))
            CALL PGRECT(REAL(K-1),REAL(K),0.,1.)
          END DO
          CALL PGSCI(1)
          CALL PGBOX('BCT',0.,0,'BC',0.,0)
        END DO
        WRITE(CDUMMY,*)YMIN0
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
          DO WHILE(CDUMMY(L:L).EQ.'0')
            L=L-1
          END DO
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGMTEXT('B',1.2,0.,0.,CDUMMY(1:L))
        END DO
        WRITE(CDUMMY,*)YMAX0
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        IF(INDEX(CDUMMY(1:L),'E').EQ.0)THEN
          DO WHILE(CDUMMY(L:L).EQ.'0')
            L=L-1
          END DO
        END IF
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGMTEXT('B',1.2,1.,1.,CDUMMY(1:L))
        END DO
C
        Y1VPORT=0.10
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C******************************************************************************
C
        SUBROUTINE PLOTSP3D(INFILE)
C Dibuja espectros en profundidad para dar sensacion de imagen tridimensional
        IMPLICIT NONE
        CHARACTER*(*) INFILE
C
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
        REAL READF
C
        INTEGER NBUT
        PARAMETER(NBUT=6)
C
        INTEGER NB
        INTEGER I,J,KK,L
        INTEGER I0
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NNC1,NNC2
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        REAL XOFFSET,YOFFSET
        REAL XMAX,XMIN,YMAX,YMAX0,YMIN,YMIN0
        REAL DX,DY
        REAL XC,YC
        REAL S(NCMAX),A(NCMAX,NSMAX)
        REAL X(NCMAX),F(NCMAX)
        REAL FMIN,FMAX
        CHARACTER*1 CH,COPC
        CHARACTER*10 BLABEL(NBUT)
        CHARACTER*40 POSICION,CDUMMY
        CHARACTER*255 LOCALGLABEL
        LOGICAL REPLOT
        LOGICAL LCOLOR(MAX_ID_RED)
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLKDEVICE1/NTERM,IDN
        COMMON/BLKDEVICE2/LCOLOR
C------------------------------------------------------------------------------
        DATA (BLABEL(KK),KK=1,NBUT)/'up','down','zoom','whole',
     +   'y-cuts','exit'/
C
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO KK=1,18
          CALL BUTTSEX(KK,.FALSE.)
        END DO
        DO KK=1,NBUT
          CALL BUTTON(KK,BLABEL(KK),0)
        END DO
C
        WRITE(CDUMMY,'(I10,A1,I10)')1,',',NCHAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'1st & last channel '
        CALL READ2I(CDUMMY(1:L),NC1,NC2)
        WRITE(CDUMMY,'(I10,A1,I10)')1,',',NSCAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        WRITE(*,100)'1st & last scan    '
        CALL READ2I(CDUMMY(1:L),NS1,NS2)
C
        YMIN0=A(NC1,NS1)
        YMAX0=YMIN0
        DO I=NS1,NS2
          DO J=NC1,NC2
            YMIN0=AMIN1(YMIN0,A(J,I))
            YMAX0=AMAX1(YMAX0,A(J,I))
          END DO
        END DO
        WRITE(*,*)
        WRITE(*,100)'YMIN = '
        WRITE(*,*)YMIN0
        WRITE(*,100)'YMAX = '
        WRITE(*,*)YMAX0
        WRITE(*,100)'Change Limits (y/n) '
        COPC(1:1)=READC('n','yn')
        IF(COPC.EQ.'y')THEN
          WRITE(*,100)'YMIN '
          WRITE(CDUMMY,*)YMIN0
          YMIN0=READF(CDUMMY)
          WRITE(*,100)'YMAX '
          WRITE(CDUMMY,*)YMAX0
          YMAX0=READF(CDUMMY)
        END IF
C
        I0=(NS1+NS2)/2
C
35      CONTINUE
C       WRITE(*,100)'X-offset between each spectrum (channels)'
C       XOFFSET=READF('@')
        XOFFSET=REAL(NC2-NC1+1)/(4*REAL(NS2-NS1+1))
C       WRITE(*,100)'Y-offset between each spectrum (counts)  '
C       YOFFSET=READF('@')
        YOFFSET=0.8*(YMAX0-YMIN0)/REAL(NS2-NS1+1)
        XMIN=REAL(NC1)-1.
        XMAX=REAL(NC2)+REAL(NS2-NS1)*XOFFSET+1.
        YMAX=YMAX0+REAL(NS2-NS1-1)*YOFFSET
        YMAX=YMAX+YOFFSET
        YMIN=YMIN0-YOFFSET
        DY=YMAX-YMIN
        DX=XMAX-XMIN
C
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          IF(ITERM.EQ.1)THEN
            CALL RPGERASW(0.,1.,0.,0.80)
ccc         CALL RPGERAS
            CALL RPGENV(XMIN,XMAX,YMIN,YMAX,0,0)
          ELSE
            CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,-2)
ccc         CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
            CALL PGBOX('BCNITS',0.0,0,'BCNITS',0.0,0)
          END IF
          CALL PGIDEN_RED
          CALL PGLABEL('channel','counts',CHAR(32))
          L=TRUELEN(INFILE)
          LOCALGLABEL(1:6)='file: '
          LOCALGLABEL(1+6:L+6)=INFILE(1:L)
          CALL PGMTEXT('T',2.0,0.5,0.5,LOCALGLABEL(1:L+6))
          IF(LCOLOR(ITERM)) CALL PGSCI(2)
        END DO
        DO I=NS2,NS1,-1
          DO J=NC1,NC2
            X(J-NC1+1)=REAL(J)+REAL(I-NS1)*XOFFSET
            S(J-NC1+1)=A(J,I)+REAL(I-NS1)*YOFFSET
          END DO
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
          END DO
        END DO
        CH=' '
        DO J=NC1,NC2
          X(J-NC1+1)=REAL(J)+REAL(I0-NS1)*XOFFSET
          S(J-NC1+1)=A(J,I0)+REAL(I0-NS1)*YOFFSET
          F(J-NC1+1)=A(J,I0)
        END DO
        CALL FINDMM(NC2-NC1+1,F,FMIN,FMAX)
        CALL PRINTD(I0,FMIN,FMAX)
        WRITE(POSICION,'(A,I4)')'Scan # ',I0
        DO ITERM=NTERM,1,-1
          CALL PGSLCT(IDN(ITERM))
          CALL PGSCI(1)
          CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
          CALL PGPTEXT(XMIN+DX/40,YMAX-DY/20,0.,0.,POSICION)
        END DO
C
        DO WHILE(CH.NE.'Q')
          CH=' '
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          CALL IFBUTTON(XC,YC,NB)
          IF(NB.EQ.0)THEN
            WRITE(*,100)'Channel: '
            WRITE(*,*)XC-REAL(I0-NS1)*XOFFSET
          ELSE
            CALL BUTTON(NB,BLABEL(NB),4)
            CALL BUTTON(NB,BLABEL(NB),1)
          END IF
          IF((NB.EQ.1).OR.(NB.EQ.2))THEN
            REPLOT=.TRUE.
          ELSE
            REPLOT=.FALSE.
          END IF
          IF(NB.EQ.1)THEN
            I0=I0+1
          END IF
          IF(NB.EQ.2)THEN
            I0=I0-1
          END IF
          IF(NB.EQ.6)THEN
            CH='Q'
          END IF
          IF(REPLOT)THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF(LCOLOR(ITERM)) CALL PGSCI(2)
              CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
              IF(LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
            IF(I0.LT.NS1) I0=NS2
            IF(I0.GT.NS2) I0=NS1
            DO J=NC1,NC2
              X(J-NC1+1)=REAL(J)+REAL(I0-NS1)*XOFFSET
              S(J-NC1+1)=A(J,I0)+REAL(I0-NS1)*YOFFSET
              F(J-NC1+1)=A(J,I0)
            END DO
            CALL FINDMM(NC2-NC1+1,F,FMIN,FMAX)
            CALL PRINTD(I0,FMIN,FMAX)
            WRITE(POSICION,'(A,I4)')'Scan # ',I0
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGBIN(NC2-NC1+1,X,S,.TRUE.)
              CALL PGSCI(0)
              CALL PGPTEXT(XMIN+DX/40,YMAX-DY/20,0.,0.,POSICION)
              CALL PGSCI(1)
              CALL PGPTEXT(XMIN+DX/40,YMAX-DY/20,0.,0.,POSICION)
            END DO
          END IF
          IF(NB.EQ.3)THEN
            WRITE(*,101)'Press region to be zoomed (first and last '//
     +       'channel)...'
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(6,0,0.,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            NNC1=INT(XC-REAL(I0-NS1)*XOFFSET+0.5)
            IF(NNC1.LT.1) NNC1=1
            IF(NNC1.GT.NCHAN) NNC1=NCHAN
            WRITE(*,110)'First channel: ',NNC1
            IF(LCOLOR(1)) CALL PGSCI(5)
            CALL RPGBAND(4,0,XC,0.,XC,YC,CH)
            IF(LCOLOR(1)) CALL PGSCI(1)
            NNC2=INT(XC-REAL(I0-NS1)*XOFFSET+0.5)    
            IF(NNC2.LT.1) NNC2=1
            IF(NNC2.GT.NCHAN) NNC2=NCHAN
            WRITE(*,110)'Last  channel: ',NNC2
            IF(NNC1.EQ.NNC2)THEN
              WRITE(*,101)'ERROR: Invalid limits.'
            ELSE
              NC1=NNC1
              NC2=NNC2
              CALL RPGERASW(0.,1.,0.,0.80)
ccc           CALL RPGERAS
              CALL BUTTON(NB,BLABEL(NB),0)
              GOTO 35
            END IF
          END IF
          IF(NB.EQ.4)THEN
            NC1=1
            NC2=NCHAN
            CALL RPGERASW(0.,1.,0.,0.80)
ccc         CALL RPGERAS
            CALL BUTTON(NB,BLABEL(NB),0)
            GOTO 35
          END IF
          IF(NB.EQ.5)THEN
            WRITE(*,100)'New YMIN '
            WRITE(CDUMMY,*)YMIN
            YMIN0=READF(CDUMMY)
            WRITE(*,100)'New YMAX '
            WRITE(CDUMMY,*)YMAX0
            YMAX0=READF(CDUMMY)
            CALL RPGERASW(0.,1.,0.,0.80)
ccc         CALL RPGERAS
            CALL BUTTON(NB,BLABEL(NB),0)
            GOTO 35
          END IF
          IF(NB.NE.0) CALL BUTTON(NB,BLABEL(NB),0)
        END DO
C
100     FORMAT(A,$)
101     FORMAT(A)
110     FORMAT(A,I6)
        END
C
C******************************************************************************
C
        SUBROUTINE PRINTD(I0,FMIN,FMAX)
        IMPLICIT NONE
        INTEGER I0
        INTEGER L,M,N
        REAL FMIN,FMAX
        CHARACTER*70 CDUMMY,CDUMMYF
C
        WRITE(CDUMMY,*)'Scan#',I0
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CDUMMYF(1:L)=CDUMMY(1:L)
        WRITE(CDUMMY,*)'Min=',FMIN
        CALL RMBLANK(CDUMMY,CDUMMY,M)
        CDUMMYF(L+1:L+5)='     '
        CDUMMYF(L+1+5:L+M+5)=CDUMMY(1:M)
        WRITE(CDUMMY,*)'Max=',FMAX
        CALL RMBLANK(CDUMMY,CDUMMY,N)
        CDUMMYF(L+M+5+1:L+M+5+5)='     '
        CDUMMYF(L+M+5+5+1:L+M+5+5+N)=CDUMMY(1:N)
        WRITE(*,'(A)')CDUMMYF(1:L+M+N+10)
        END
C
C******************************************************************************
C Calcula (usando las imágenes auxiliares) la diferencia entre el pixel
C considerado en la imagen principal y el valor promedio en el mismo pixel
C para las imagenes auxiliares.
        SUBROUTINE CHECK_AUXFRAME(J,I,DEV)
        IMPLICIT NONE
        INTEGER J,I
        REAL DEV
C
        INCLUDE 'redlib.inc'
C
        INTEGER NAUXMAX
        PARAMETER (NAUXMAX=3)
        INTEGER NPPMAX
        PARAMETER (NPPMAX=101)
C
        REAL FMEDIAN1
C
        INTEGER K,II,JJ,II_,JJ_
        INTEGER NAUX
        INTEGER NAUXFRAME
        INTEGER NC1,NC2,NS1,NS2
        INTEGER NFIT,NPP,NDP
        REAL A(NCMAX,NSMAX)
        REAL AUXFRAME(NCMAX,NSMAX,NAUXMAX)
        REAL AREPLACE(NPPMAX,NPPMAX)
        REAL FMEDIAN(0:NAUXMAX)
        REAL PIXEL(NPPMAX*NPPMAX)
        REAL DIFF
C
        COMMON/BLKDATA/A,NSCAN,NCHAN
        COMMON/BLK5/NFIT,NPP,NDP
        COMMON/BLK17A/NAUXFRAME
        COMMON/BLK17B/AUXFRAME
        COMMON/BLK23/AREPLACE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C Protecciones
        IF(NAUXFRAME.LT.1)THEN
          WRITE(*,100) 'Number of auxiliary frames: '
          WRITE(*,*) NAUXFRAME
          WRITE(*,101) 'ERROR: LOTHERCR cannot be computed.'
          RETURN
        END IF
C------------------------------------------------------------------------------
C Coordenadas de las esquinas de una ventana de lado NPP centrada en
C el pixel (J,I)
        NC1=J-NPP/2
        IF(NC1.LT.1) NC1=1
        NC2=NC1+(NPP-1)
        IF(NC2.GT.NCHAN)THEN
          NC2=NCHAN
          NC1=NCHAN-(NPP-1)
        END IF
        NS1=I-NPP/2
        IF(NS1.LT.1) NS1=1
        NS2=NS1+(NPP-1)
        IF(NS2.GT.NSCAN)THEN
          NS2=NSCAN
          NS1=NSCAN-(NPP-1)
        END IF
C Calculamos la mediana en cada imagen auxiliar dentro de la ventana de 
C lado NPP pixels
        DO NAUX=1,NAUXFRAME
          K=0
          DO II=NS1,NS2
            DO JJ=NC1,NC2
              K=K+1
              PIXEL(K)=AUXFRAME(JJ,II,NAUX)
            END DO
          END DO
          FMEDIAN(NAUX)=FMEDIAN1(K,PIXEL)
        END DO
C Calculamos la mediana en la imagen principal y almacenamos el valor en
C la variable FMEDIAN(0).
        K=0
        DO II=NS1,NS2
          DO JJ=NC1,NC2
            K=K+1
            PIXEL(K)=A(JJ,II)
          END DO
        END DO
        FMEDIAN(0)=FMEDIAN1(K,PIXEL)
C Generamos una subimagen (NPP*NPP) con la media de las imagenes auxiliares
C (escaladas antes a la señal de la imagen principal)
        DO II=NS1,NS2
          II_=II-NS1+1
          DO JJ=NC1,NC2
            JJ_=JJ-NC1+1
            AREPLACE(JJ_,II_)=0.
            DO NAUX=1,NAUXFRAME
              AREPLACE(JJ_,II_)=AREPLACE(JJ_,II_)+
     +         AUXFRAME(JJ,II,NAUX)*FMEDIAN(0)/FMEDIAN(NAUX)
            END DO
            AREPLACE(JJ_,II_)=AREPLACE(JJ_,II_)/REAL(NAUXFRAME)
          END DO
        END DO
C Calculamos la mediana de las diferencias entre la imagen principal y el
C promedio anterior. Esto nos sirve como estimacion de la dispersión entre 
C ambas imagenes.
        K=0
        DO II=NS1,NS2
          II_=II-NS1+1
          DO JJ=NC1,NC2
            JJ_=JJ-NC1+1
            K=K+1
            PIXEL(K)=ABS(A(JJ,II)-AREPLACE(JJ_,II_))
          END DO
        END DO
        DIFF=FMEDIAN1(K,PIXEL)
C
        IF(DIFF.EQ.0.0)THEN
          DEV=0.0
        ELSE
          DEV=(A(J,I)-AREPLACE(J-NC1+1,I-NS1+1))/DIFF
        END IF
C
100     FORMAT(A,$)
101     FORMAT(A)
        END
