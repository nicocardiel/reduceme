!------------------------------------------------------------------------------
! Version 24-November-2026                                  file: lkeepcard.f90
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
! LOGICAL FUNCTION LKEEPCARD(CARD)
!
! Input: CARD
! Output: LKEEPCARD
!
! Returns .TRUE. if the FITS header card CARD of an input FITS file must be
! copied into the FITS output files. The following cards are not copied:
! - structural keywords (SIMPLE, BITPIX, NAXISn, EXTEND, BSCALE, BZERO,
!   BLANK, XTENSION, EXTNAME, CHECKSUM, END...), written again by CFITSIO;
! - keywords written by REDUCEME from its own header (OBJECT, EXPTIME,
!   AIRMASS, ERRFRAME, FITSFILE, REDCOMM, DATE, LONGSTRN) and the old ones
!   of writefits (TIMEXPO, FITSFI);
! - all the world coordinate system keywords (CTYPEn, CRPIXn, CRVALn,
!   CDELTn, CDi_j, PCi_j, PVi_m, WCSAXES..., including alternate
!   descriptions with a final letter, and the IRAF ones: LTVn, LTMi_j,
!   WATi_nnn, DC-FLAG, DISPAXIS), which may not be valid after the
!   processing (REDUCEME writes those of the wavelength axis);
! - the COMMENT cards that CFITSIO adds to every new file (FITS definition
!   and long string convention).
! CONTINUE cards must be handled by the caller (they belong to the previous
! card).
!
! CHARACTER*80 CARD -> FITS header card
!
!Comment
!------------------------------------------------------------------------------
        LOGICAL FUNCTION LKEEPCARD(CARD)
        IMPLICIT NONE
        CHARACTER*(*) CARD
!
        INTEGER NEXACT
        PARAMETER (NEXACT=31)
        INTEGER NINDEX
        PARAMETER (NINDEX=11)
        INTEGER NIJ
        PARAMETER (NIJ=6)
        INTEGER I
        CHARACTER*8 KEY
        CHARACTER*8 CEXACT(NEXACT)
        CHARACTER*8 CINDEX(NINDEX)
        CHARACTER*8 CIJ(NIJ)
        DATA CEXACT/'SIMPLE  ','BITPIX  ','NAXIS   ','EXTEND  ',&
     &   'BSCALE  ','BZERO   ','BLANK   ','PCOUNT  ','GCOUNT  ',&
     &   'XTENSION','EXTNAME ','EXTVER  ','EXTLEVEL','INHERIT ',&
     &   'CHECKSUM','DATASUM ','END     ','DATE    ','LONGSTRN',&
     &   'OBJECT  ','EXPTIME ','TIMEXPO ','AIRMASS ','ERRFRAME',&
     &   'FITSFILE','FITSFI  ','REDCOMM ','LONPOLE ','LATPOLE ',&
     &   'DC-FLAG ','DISPAXIS'/
! keywords followed by an axis number (and an optional alternate letter)
        DATA CINDEX/'NAXIS','CTYPE','CUNIT','CRPIX','CRVAL','CDELT',&
     &   'CROTA','CRDER','CSYER','CNAME','LTV'/
! keywords followed by two indices i_j (and an optional alternate letter)
        DATA CIJ/'CD','PC','PV','PS','LTM','WAT'/
!------------------------------------------------------------------------------
        KEY=CARD(1:8)
        LKEEPCARD=.FALSE.
!
        DO I=1,NEXACT
          IF(KEY.EQ.CEXACT(I)) RETURN
        END DO
        IF(KEY(1:7).EQ.'WCSAXES') RETURN
        IF(KEY(1:7).EQ.'WCSNAME') RETURN
        DO I=1,NINDEX
          IF(LINDEXED(KEY,CINDEX(I)(1:TRUELEN(CINDEX(I))),.FALSE.)) RETURN
        END DO
        DO I=1,NIJ
          IF(LINDEXED(KEY,CIJ(I)(1:TRUELEN(CIJ(I))),.TRUE.)) RETURN
        END DO
! COMMENT cards added by CFITSIO to every new file
        IF(KEY.EQ.'COMMENT ')THEN
          IF(INDEX(CARD,'FITS (Flexible Image Transport System) format').GT.0)&
     &     RETURN
          IF(INDEX(CARD,'and Astrophysics'', volume 376, page 359').GT.0)&
     &     RETURN
          IF(INDEX(CARD,'This FITS file may contain long string').GT.0)&
     &     RETURN
          IF(INDEX(CARD,'continued over multiple keywords.').GT.0) RETURN
          IF(INDEX(CARD,'character at the end of each substring').GT.0)&
     &     RETURN
          IF(INDEX(CARD,'on the next keyword which has the name').GT.0)&
     &     RETURN
        END IF
        LKEEPCARD=.TRUE.
!
        CONTAINS
!
! .TRUE. if KEY is CPREFIX followed by a number (LIJ=.FALSE.) or by two
! numbers separated by '_' (LIJ=.TRUE.), with an optional final letter
! (alternate world coordinate system descriptions)
        LOGICAL FUNCTION LINDEXED(KEY,CPREFIX,LIJ)
        IMPLICIT NONE
        CHARACTER*8 KEY
        CHARACTER*(*) CPREFIX
        LOGICAL LIJ
!
        INTEGER LP,LK,J,NUNDER,NDIGITS
        CHARACTER*1 C
!
        LINDEXED=.FALSE.
        LP=LEN(CPREFIX)
        LK=TRUELEN(KEY)
        IF(LK.LE.LP) RETURN
        IF(KEY(1:LP).NE.CPREFIX) RETURN
! optional final alternate letter
        C=KEY(LK:LK)
        IF((C.GE.'A').AND.(C.LE.'Z')) LK=LK-1
        IF(LK.LE.LP) RETURN
        NUNDER=0
        NDIGITS=0
        DO J=LP+1,LK
          C=KEY(J:J)
          IF((C.GE.'0').AND.(C.LE.'9'))THEN
            NDIGITS=NDIGITS+1
          ELSEIF(C.EQ.'_')THEN
            NUNDER=NUNDER+1
            IF((J.EQ.LP+1).OR.(J.EQ.LK)) RETURN
          ELSE
            RETURN
          END IF
        END DO
        IF(NDIGITS.EQ.0) RETURN
        IF(LIJ)THEN
          LINDEXED=(NUNDER.EQ.1)
        ELSE
          LINDEXED=(NUNDER.EQ.0)
        END IF
        END FUNCTION LINDEXED
!
        END
