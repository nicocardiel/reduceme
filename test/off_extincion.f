C Simula el efecto de errores en la aplicacion de la curva de extincion
C en la observacion del D4000
	PROGRAM EXTINCION
	IMPLICIT NONE
	INCLUDE 'futils.inc'
C
	INTEGER NMAX
	PARAMETER (NMAX=20)
C
	INTEGER I
	INTEGER NSEED
	REAL FA_REAL(NMAX),FR_REAL(NMAX),D4000_REAL(NMAX)
	REAL FA_EXTI(NMAX),FR_EXTI(NMAX),D4000_EXTI(NMAX)
	REAL FA_CORR(NMAX),FR_CORR(NMAX),D4000_CORR(NMAX)
	REAL AIRMASS(NMAX),ZENIT,ZENITMAX
	REAL AMAGA,AMAGR,ERRAMAGA,ERRAMAGR,ERRAMAGAF,ERRAMAGRF
	REAL ERRSYSAMAGA,ERRSYSAMAGR
	REAL XMIN,XMAX
	REAL RANRED
C------------------------------------------------------------------------------
	CALL PGBEGIN(0,'?',1,1)
	CALL PGSCH(1.5)
	CALL PGASK(.FALSE.)
C------------------------------------------------------------------------------
	AMAGA=0.3440              !extincion (mag/airmass) en 3850 A (La Palma)
	AMAGR=0.2519              !extincion (mag/airmass) en 4150 A (La Palma)
	NSEED=-1
	WRITE(*,100) 'Maximum zenital distance (degrees) '
	ZENITMAX=READF('60.0')
	WRITE(*,101) 'EXTINCTION CURVE: '
	WRITE(*,100) 'BLUE: maximum random error (mag) '
	ERRAMAGA=READF('0.10')
	WRITE(*,100) 'BLUE: systematic error (mag) '
	ERRSYSAMAGA=READF('0.05')
	WRITE(*,100) 'RED: maximum random error (mag) '
	ERRAMAGR=READF('0.05')
	WRITE(*,100) 'RED: systematic error (mag) '
	ERRSYSAMAGR=READF('0.00')
	DO I=1,NMAX
	  FA_REAL(I)=1.0
	  ZENIT=ZENITMAX*RANRED(NSEED)    !distancia zenital entre 0 y ZENITMAX
	  AIRMASS(I)=1./COS(ZENIT*3.141593/180.)
	  FR_REAL(I)=FA_REAL(I)+2.5*RANRED(NSEED) !creamos B4000 entre 1. y 3.5
	  D4000_REAL(I)=1.1619*FR_REAL(I)/FA_REAL(I)
	  FA_EXTI(I)=FA_REAL(I)/10.**(0.4*AIRMASS(I)*AMAGA)
	  FR_EXTI(I)=FR_REAL(I)/10.**(0.4*AIRMASS(I)*AMAGR)
	  D4000_EXTI(I)=1.1619*FR_EXTI(I)/FA_EXTI(I)
	  ERRAMAGAF=ERRSYSAMAGA+AIRMASS(I)*ERRAMAGA*RANRED(NSEED)
	  ERRAMAGRF=ERRSYSAMAGR+AIRMASS(I)*ERRAMAGR*RANRED(NSEED)
	  FA_CORR(I)=FA_EXTI(I)*10.**(0.4*AIRMASS(I)*(AMAGA+ERRAMAGAF))
	  FR_CORR(I)=FR_EXTI(I)*10.**(0.4*AIRMASS(I)*(AMAGR+ERRAMAGRF))
	  D4000_CORR(I)=1.1619*FR_CORR(I)/FA_CORR(I)
	END DO
	XMIN=0.8
	XMAX=4.8
	CALL PGENV(XMIN,XMAX,XMIN,XMAX,1,0)
	CALL PGLABEL('D\\d4000\\u real','D\\d4000\\u with extinction',
     +   ' ')
	CALL PGPOINT(NMAX,D4000_REAL,D4000_EXTI,17)
	CALL PGSCI(3)
	CALL PGPOINT(NMAX,D4000_REAL,D4000_CORR,5)
	CALL PGSCI(2)
	CALL PGMOVE(XMIN,XMIN)
	CALL PGDRAW(XMAX,XMAX)
	CALL PGSCI(1)
C------------------------------------------------------------------------------
	CALL PGEND
C------------------------------------------------------------------------------
	STOP
100	FORMAT(A,$)
101	FORMAT(A)
	END
