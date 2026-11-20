!------------------------------------------------------------------------------
! Version 20-November-2026                                   File: button_m.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This file is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
! Module BUTTON_M: graphic buttons (PGPLOT).
! The procedures are kept in their own files, included below; each one
! declares IMPLICIT NONE.
!------------------------------------------------------------------------------
        MODULE BUTTON_M
        USE FUTILS_M
! global variables of the graphic buttons (previously in button.inc,
! in the COMMON blocks /BUTTON1/ to /BUTTON8/)
        INTEGER MAX_NBUTT
        PARAMETER (MAX_NBUTT=1000)
        INTEGER MAX_ID_BUTT
        PARAMETER (MAX_ID_BUTT=8)
        INTEGER MAX_XBUTT,MAX_YBUTT
        INTEGER PGSCF_BUTT
        REAL PGSCH_BUTT
        REAL YTEXT_BUTT
        REAL X1VPORT,X2VPORT
        REAL X3VPORT,X4VPORT
        REAL Y1VPORT,Y2VPORT
        REAL Y3VPORT,Y4VPORT
        LOGICAL EXIST_BUTT(MAX_NBUTT)
        LOGICAL ITICKS_BUTT
        LOGICAL MODOTEXT_BUTT
        LOGICAL MODOTEXT_PLOTBUTT
        LOGICAL MODOTEXT_HTML
! these variables are only used by the procedures of this module
        PRIVATE :: MAX_NBUTT,MAX_ID_BUTT,MAX_XBUTT,MAX_YBUTT,PGSCF_BUTT,PGSCH_BUTT,&
     &   YTEXT_BUTT,X1VPORT,X2VPORT,X3VPORT,X4VPORT,&
     &   Y1VPORT,Y2VPORT,Y3VPORT,Y4VPORT,EXIST_BUTT,&
     &   ITICKS_BUTT,MODOTEXT_BUTT,MODOTEXT_PLOTBUTT,MODOTEXT_HTML
        CONTAINS
        INCLUDE 'buttmis.f90'
        INCLUDE 'button.f90'
        INCLUDE 'buttqbr.f90'
        INCLUDE 'buttqcf.f90'
        INCLUDE 'buttqch.f90'
        INCLUDE 'buttqex.f90'
        INCLUDE 'buttqit.f90'
        INCLUDE 'buttqpr.f90'
        INCLUDE 'buttqxb.f90'
        INCLUDE 'buttqyb.f90'
        INCLUDE 'buttqytext.f90'
        INCLUDE 'buttsbr.f90'
        INCLUDE 'buttscf.f90'
        INCLUDE 'buttsch.f90'
        INCLUDE 'buttsex.f90'
        INCLUDE 'buttsit.f90'
        INCLUDE 'buttspr.f90'
        INCLUDE 'buttsxb.f90'
        INCLUDE 'buttsyb.f90'
        INCLUDE 'buttsytext.f90'
        INCLUDE 'ifbutton.f90'
        INCLUDE 'rpgband.f90'
        INCLUDE 'rpgbegin.f90'
        INCLUDE 'rpgbegok.f90'
        INCLUDE 'rpgenv.f90'
        INCLUDE 'rpgeras.f90'
        INCLUDE 'rpgerasb.f90'
        INCLUDE 'rpgerasw.f90'
        END MODULE BUTTON_M
