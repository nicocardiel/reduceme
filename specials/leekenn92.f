	PROGRAM LEEKENN92
	IMPLICIT NONE
C
	INTEGER I,J,J1,J2,K
	INTEGER NCHAR
	INTEGER NSCAN,NCHAN
	INTEGER TRUELEN
	REAL STWV,DISP
	REAL AIRMASS,TIMEXPOS
	REAL SPECTRA(1726,56)
	CHARACTER*12 IDENTIFICATION
	CHARACTER*255 OBJECT
	CHARACTER*255 FITSFILE
	CHARACTER*255 COMMENT
C------------------------------------------------------------------------------
C read input file
	OPEN(10,FILE='catalog.dat',STATUS='OLD',FORM='FORMATTED')
	DO I=1,56
	  READ(10,*) !skip first line
	  DO K=1,172
	    J1=(K-1)*10+1
	    J2=J1+9
	    READ(10,*) (SPECTRA(J,I),J=J1,J2)
	  END DO
	  J1=1721
	  J2=1726
	  READ(10,*) (SPECTRA(J,I),J=J1,J2)
	END DO
	CLOSE(10)
C------------------------------------------------------------------------------
C write output file
	OPEN(10,FILE='catalog.u',STATUS='NEW',FORM='UNFORMATTED')
	IDENTIFICATION='abcdefghijkl'
	WRITE(10) IDENTIFICATION
	NSCAN=56
	NCHAN=1726
	WRITE(10) NSCAN,NCHAN
	STWV=3650.0
	DISP=2.0
	WRITE(10) STWV,DISP
	AIRMASS=0.
	WRITE(10) AIRMASS
	TIMEXPOS=0.
	WRITE(10) TIMEXPOS
	OBJECT='KENNICUTT 1992'
	NCHAR=TRUELEN(OBJECT)
	WRITE(10) NCHAR
	IF(NCHAR.GT.0) WRITE(10) OBJECT(1:NCHAR)
	FITSFILE='NONE'
	NCHAR=TRUELEN(FITSFILE)
	WRITE(10) NCHAR
	IF(NCHAR.GT.0) WRITE(10) FITSFILE(1:NCHAR)
	COMMENT='Spectrophotometric Atlas of Galaxies'
	NCHAR=TRUELEN(COMMENT)
	WRITE(10) NCHAR
	IF(NCHAR.GT.0) WRITE(10) COMMENT(1:NCHAR)
	DO I=1,NSCAN
	  WRITE(10) (SPECTRA(J,I),J=1,NCHAN)
	END DO
	CLOSE(10)
C------------------------------------------------------------------------------
C end of program
	STOP
	END
C
C******************************************************************************
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
