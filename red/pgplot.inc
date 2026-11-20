!-----------------------------------------------------------------------
! PGPLOT: common block definition.
!-----------------------------------------------------------------------
! Maximum number of concurrent devices (should match GRIMAX).
!-----------------------------------------------------------------------
      INTEGER PGMAXD
      PARAMETER (PGMAXD=8)
!-----------------------------------------------------------------------
! Indentifier of currently selected device.
!-----------------------------------------------------------------------
      INTEGER PGID
!-----------------------------------------------------------------------
! Device status (indexed by device identifier).
!-----------------------------------------------------------------------
! PGDEVS  =0 if device is not open; 1 if device is open.
! PGADVS  Set to 0 by PGBEGIN, set to 1 by PGPAGE; used to suppress
!         the prompt for the first page.
! PROMPT  If .TRUE., ask user before clearing page; set by PGASK
!         and (indirectly) by PGBEGIN, used in PGENV.
! PGBLEV  Buffering level: incremented by PGBBUF, decremented by
!         PGEBUF.
! PGPFIX  TRUE if PGPAP has been called, FALSE otherwise.
!
      INTEGER PGDEVS(PGMAXD), PGADVS(PGMAXD), PGBLEV(PGMAXD)
      LOGICAL PGPRMP(PGMAXD), PGPFIX(PGMAXD)
!-----------------------------------------------------------------------
! Panel parameters (indexed by device identification).
!-----------------------------------------------------------------------
! NX      Number of panels in x direction
! NY      Number of panels in y direction
! NXC     Ordinal number of current X panel
! NYC     Ordinal number of current Y panel
! XSZ     X dimension of panel (device units)
! YSZ     Y dimension of panel (device units)
! PGROWS  TRUE if panels are used in row order, FALSE for column
!         order.
!
      INTEGER PGNX  (PGMAXD), PGNY  (PGMAXD)
      INTEGER PGNXC (PGMAXD), PGNYC (PGMAXD)
      REAL    PGXSZ (PGMAXD), PGYSZ (PGMAXD)
      LOGICAL PGROWS(PGMAXD)
!-----------------------------------------------------------------------
! Attributes (indexed by device identification).
!-----------------------------------------------------------------------
! PGCLP   clipping enabled/disabed
! PGFAS   fill-area style
! PGCHSZ  character height
! PGAHS   arrow-head fill style
! PGAHA   arrow-head angle
! PGAHV   arrow-head vent
! PGTBCI  text background color index
! PGMNCI  lower range of color indices available to PGGRAY/PGIMAG
! PGMXCI  upper range of color indices available to PGGRAY/PGIMAG
! PGITF   type of transfer function used by PGGRAY/PGIMAG
! PGHSA   hatching line angle
! PGHSS   hatching line separation
! PGHSP   hatching line phase
!
      INTEGER PGCLP (PGMAXD)
      INTEGER PGFAS (PGMAXD)
      REAL    PGCHSZ(PGMAXD)
      INTEGER PGAHS (PGMAXD)
      REAL    PGAHA (PGMAXD)
      REAL    PGAHV (PGMAXD)
      INTEGER PGTBCI(PGMAXD)
      INTEGER PGMNCI(PGMAXD)
      INTEGER PGMXCI(PGMAXD)
      INTEGER PGITF (PGMAXD)
      REAL    PGHSA (PGMAXD)
      REAL    PGHSS (PGMAXD)
      REAL    PGHSP (PGMAXD)
!-----------------------------------------------------------------------
! Viewport parameters (indexed by device identification); all are device
! coordinates:
!-----------------------------------------------------------------------
! PGXOFF  X coordinate of blc of viewport.
! PGYOFF  Y coordinate of blc of viewport.
! PGXVP   X coordinate of blc of viewport, relative to blc of subpage.
! PGYVP   Y coordinate of blc of viewport, relative to blc of subpage.
! PGXLEN  Width of viewport. 
! PGYLEN  Height of viewport.
!
      REAL   PGXOFF(PGMAXD), PGYOFF(PGMAXD)
      REAL   PGXVP (PGMAXD), PGYVP (PGMAXD)
      REAL   PGXLEN(PGMAXD), PGYLEN(PGMAXD)
!-----------------------------------------------------------------------
! Scaling parameters (indexed by device identification):
!-----------------------------------------------------------------------
! PGXORG  device coordinate value corresponding to world X=0
! PGYORG  device coordinate value corresponding to world Y=0
! PGXSCL  scale in x (device units per world coordinate unit)
! PGYSCL  scale in y (device units per world coordinate unit)
! PGXPIN  device x scale in device units/inch
! PGYPIN  device y scale in device units/inch
! PGXSP   Character X spacing (device units)
! PGYSP   Character Y spacing (device units)
!
      REAL   PGXORG(PGMAXD), PGYORG(PGMAXD)
      REAL   PGXSCL(PGMAXD), PGYSCL(PGMAXD)
      REAL   PGXPIN(PGMAXD), PGYPIN(PGMAXD)
      REAL   PGXSP (PGMAXD), PGYSP (PGMAXD)
!-----------------------------------------------------------------------
! Window parameters (indexed by device identification); all are world
! coordinate values:
!-----------------------------------------------------------------------
! PGXBLC  world X at bottom left corner of window
! PGXTRC  world X at top right corner of window
! PGYBLC  world Y at bottom left corner of window
! PGYTRC  world Y at top right corner of window
!
      REAL    PGXBLC(PGMAXD), PGXTRC(PGMAXD)
      REAL    PGYBLC(PGMAXD), PGYTRC(PGMAXD)
!-----------------------------------------------------------------------
! The following parameters are used in the contouring routines to pass
! information to the action routine. They do not need to be indexed.
!-----------------------------------------------------------------------
! TRANS   Transformation matrix for contour plots; copied
!         from argument list by PGCONT and used by PGCP.
!
      INTEGER PGCINT, PGCMIN
      REAL    TRANS(6)
      CHARACTER*32 PGCLAB
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      COMMON  /PGPLT1/ PGID,PGDEVS,PGADVS,PGNX,  PGNY,  PGNXC, PGNYC ,&
     &        PGXPIN,PGYPIN,PGXSP, PGYSP, PGXSZ, PGYSZ,&
     &        PGXOFF,PGYOFF,PGXVP, PGYVP, PGXLEN,PGYLEN,PGXORG,PGYORG,&
     &        PGXSCL,PGYSCL,PGXBLC,PGXTRC,PGYBLC,PGYTRC,TRANS,&
     &        PGPRMP,PGCLP, PGFAS, PGCHSZ,PGBLEV,PGROWS,&
     &        PGAHS, PGAHA, PGAHV, PGTBCI,PGMNCI,PGMXCI,PGCINT,PGCMIN,&
     &        PGPFIX,PGITF, PGHSA, PGHSS, PGHSP
      COMMON  /PGPLT2/ PGCLAB
      SAVE    /PGPLT1/
      SAVE    /PGPLT2/
!-----------------------------------------------------------------------
