        PROGRAM SIMPLE_WRITE
        IMPLICIT NONE
C
        INTEGER NCHAR
        INTEGER NSCAN,NCHAN
        INTEGER TRUELEN
        REAL STWV,DISP
        REAL AIRMASS,TIMEXPOS
        REAL SPECTRA(1124,1124)
        CHARACTER*12 IDENTIFICATION
        CHARACTER*255 OBJECT
        CHARACTER*255 FITSFILE
        CHARACTER*255 COMMENT
C
C open file
        OPEN(10,FILE='file000.dat',STATUS='NEW',FORM='UNFORMATTED')
C write header information
        IDENTIFICATION='abcdefghijkl'
        WRITE(10) IDENTIFICATION
        WRITE(10) NSCAN,NCHAN
        WRITE(10) STWV,DISP
        WRITE(10) AIRMASS
        WRITE(10) TIMEXPOS
        NCHAR=TRUELEN(OBJECT)
        WRITE(10) NCHAR
        IF(NCHAR.GT.0) WRITE(10) OBJECT(1:NCHAR)
        NCHAR=TRUELEN(FITSFILE)
        WRITE(10) NCHAR
        IF(NCHAR.GT.0) WRITE(10) FITSFILE(1:NCHAR)
        NCHAR=TRUELEN(COMMENT)
        WRITE(10) NCHAR
        IF(NCHAR.GT.0) WRITE(10) COMMENT(1:NCHAR)
C write data frame
        DO I=1,NSCAN
          WRITE(10) (SPECTRA(J,I),J=1,NCHAN)
        END DO
        CLOSE(10)
C end of program
        STOP
        END
C------------------------------------------------------------------------------
C return the "true" LEN of CADENA (ignoring control and blank characters)
        INTEGER FUNCTION TRUELEN(CADENA)
        IMPLICIT NONE
        CHARACTER*(*) CADENA
        INTEGER I,L
C
        L=LEN(CADENA)
        DO I=L,1,-1
          IF(ICHAR(CADENA(I:I)).GT.32)THEN
            TRUELEN=I
            RETURN
          END IF
        END DO
        TRUELEN=0
        END
