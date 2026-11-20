!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: pgiden_red.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE PGIDEN_RED
!
! Input (COMMON): THISPROGRAM,CREDUCEVERSION
!
! Write the current program name, the username, date, time and REDUCEME
! version at the bottom of the plot.
!
!Comment
!------------------------------------------------------------------------------
      SUBROUTINE PGIDEN_RED
      IMPLICIT NONE
      INTEGER TRUELEN
      INCLUDE 'redlib.inc'
!
      INCLUDE 'pgplot.inc'
      INTEGER L, M, CF, CI, LW
      INTEGER LL,LV,LPGPVER
      CHARACTER*64 TEXT,LOGNAME
      CHARACTER*255 ETIQUETAFIN,PGPVER
      REAL D, CH
!------------------------------------------------------------------------------
      CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
      CALL PGBBUF
!
! Get information for annotation.
!
      CALL GET_ENVIRONMENT_VARIABLE('LOGNAME',LOGNAME) !si LOGNAME esta definido, es mejor
      LL=TRUELEN(LOGNAME)
      IF(LL.EQ.0)THEN   !de lo contrario usamos GRUSER (aunque falla a veces!!)
        CALL GRUSER(TEXT, L)
      ELSE
        L=LL
        TEXT(1:L)=LOGNAME(1:L)
      END IF
      TEXT(L+1:) = ' '
!
      CALL GRDATE(TEXT(L+2:), M)
      L = L+1+M
!
! Save current attributes.
!
      CALL PGQCF(CF)
      CALL PGQCI(CI)
      CALL PGQLW(LW)
      CALL PGQCH(CH)
!
! Change attributes and write text.
!
      CALL PGSCF(3)
      CALL PGSCI(1)
      CALL PGSLW(1)
      CALL PGSCH(0.6*REAL(PGNY(PGID)))
      LL=TRUELEN(THISPROGRAM)
      ETIQUETAFIN(1:LL)=THISPROGRAM(1:LL)
      LL=LL+1
      ETIQUETAFIN(LL:LL+1)=' ('
      LL=LL+2
      LV=LEN(CREDUCEVERSION)
      ETIQUETAFIN(LL:LL+LV-1)=CREDUCEVERSION
      LL=LL+LV
      ETIQUETAFIN(LL:LL+8)=')  user:'
      LL=LL+9
      ETIQUETAFIN(LL:LL+L-1)=TEXT(1:L)
      LL=LL+L-1
      CALL GRLEN(ETIQUETAFIN(1:LL),D)
      CALL GRTEXT(.FALSE., 0.0, .TRUE.,&
     &            PGXSZ(PGID)-D-2.0+PGXSZ(PGID)*REAL(PGNX(PGID)-1),&
     &            2.0+PGYSZ(PGID)/130.0,&
     &            ETIQUETAFIN(1:LL))
      ETIQUETAFIN(1:12)='Using PGPLOT'
      CALL PGQINF('VERSION',PGPVER,LPGPVER)
      ETIQUETAFIN(13:13+LPGPVER-1)=PGPVER(1:LPGPVER)
      LL=13+LPGPVER-1
      CALL GRLEN(ETIQUETAFIN(1:LL),D)
      CALL GRTEXT(.FALSE., 0.0, .TRUE.,& 
     &            2.0,&
     &            2.0+PGYSZ(PGID)/130.0,& 
     &            ETIQUETAFIN(1:LL))
!
! Restore attributes.
!
      CALL PGSCF(CF)
      CALL PGSCI(CI)
      CALL PGSLW(LW)
      CALL PGSCH(CH)
      CALL PGEBUF
!
      END
