C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE  FITL(X,Y,NX,sig,iw,xmn,xmx,A,B,SIGA,SIGB,STD)
C
C------------------------------------------------------------------------C
C Linear Fit Routines                          J. Jesus Gonzalez G.
C------------------------------------------------------------------------C
C------------------------------------------------------------------------C
C                      Fits the line y = b*x + a
C     THIS ROUTINE ITERATES ELIMINATING HIGHLY DEVIANT POINTS
C     INPUT:  X,   Y - Data arrays
C                SIG - Y-error of points (<=0 if a point is to be rejected)
C                 IW - =0 if unweighted fit, weighted fit otherwise.
C            xmn,xmx - Limits of fit.
C
C     OUTPUT: 
C            B,    A - Slope and zero-ordinate.
C         SIGB, SIGA - Stimated errors on B and A.
C                STD - Unbiased Std-Deviation from the fit.
C
C     Remeber how to compute errors of predicted values:
C     VAR(y(x)) = VAR(a) + VAR(b)*(x^2 - 2*xm*x), since the
C     ab-error covariance is COV(ab)=-xm*VAR(b)
C
C------------------------------------------------------------------------C
Comment
C------------------------------------------------------------------------------
      SUBROUTINE FITL(X,Y,NX,sig,iw,xmn,xmx,A,B,SIGA,SIGB,STD)
      DIMENSION X(NX),Y(NX),sig(NX)

      NITER = 2
      ITER = 0
      NTD = 0

!10    sw=0.
      sw=0.
      SX=0.
      SY=0.
      ST2=0.
      B=0.
      IF (iw.NE.0.0) THEN
        DO 11 I=1,NX
         if (sig(I).le.0. .or.x(i).lt.xmn.or.x(i).gt.xmx) goto 11
          wt = (1.0/sig(i)**2)
          sw=sw+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
        XM=SX/sw
        DO 12 I=1,NX
          if (sig(i).le.0. .or.x(i).lt.xmn.or.x(i).gt.xmx) goto 12
          Tx=(X(I)-XM)/sig(i)
          ST2=ST2+Tx*Tx
          B=B+Tx*Y(I)/sig(i)
12      CONTINUE
      else
        DO 13 I=1,NX
          if (x(i).lt.xmn .or. x(i).gt.xmx) goto 13
          sw = sw + 1.0
          SX=SX+X(I)
          SY=SY+Y(I)
13      CONTINUE
        XM=SX/sw
        DO 14 I=1,NX
          if (x(i).lt.xmn .or. x(i).gt.xmx) goto 14
          Tx=(X(I)-XM)
          ST2=ST2+Tx*Tx
          B=B+Tx*Y(I)
14      CONTINUE
      end if
      B=B/ST2
      A=(SY-SX*B)/sw
      SIGA=SQRT((1.+SX*SX/(sw*ST2))/sw)
      SIGB=SQRT(1./ST2)

C--   -- Estimate the Unbiased Standard-deviation.
      std = 0.0
      IF (iw.NE.0.0) THEN
        sw2 = 0.0
        DO I=1,NX
          if (sig(i).gt.0. .and. x(i).ge.xmn.and.x(i).le.xmx) then
            sw2 = sw2 + (1.0/sig(i)**4)
            std = std + ((y(i)-a-b*x(i))/sig(i))**2
          end if
        end do
        pts = sw**2/sw2
        std = sqrt(std/sw)
        if (pts.gt.2.5) std = sqrt(pts/(pts-2.0))*std
      ELSE
        do I=1,NX
         if (x(i).ge.xmn.and.x(i).le.xmx) std=std+(y(i)-a-b*x(i))**2
        end do
        pts = sw
        std = sqrt(std/(pts-2.0))
      END IF

C--   -- Scale errors (assuming a good fit) if observed > expected error.
      if (pts.le.3.) then
       scale = 1.0
      else
       scale = std*sqrt(sw/pts)
      end if
      siga = siga*scale
      sigb = sigb*scale
      RETURN
      END

      SUBROUTINE FITLT(X,Y,NX,sig,iw,xmn,xmx,A,B,SIGA,SIGB,
     &   std,XM,T,ntd)
C------------------------------------------------------------------------C
C                      Fits the line y = b*x + a
C     THIS ROUTINE ITERATES ELIMINATING HIGHLY DEVIANT POINTS
C     INPUT:  X,   Y - Data arrays
C                SIG - Y-error of points (<=0 if a point is to be rejected)
C                 IW - =0 if unweighted fit, weighted fit otherwise.
C            xmn,xmx - Limits of fit.
C                  T - Rejection threshold (t-StdDeviations or more).
C
C     OUTPUT: 
C            B,    A - Slope and zero-ordinate.
C         SIGB, SIGA - Stimated errors on B and A.
C                STD - Unbiased Std-Deviation from the fit.
C                 XM - Mean abscisa
C                NTD - Total Number of reject points. 
C             SIG(i) - returns negative if i-th point was rejected.
C
C     Remeber how to compute errors of predicted values:
C     VAR(y(x)) = VAR(a) + VAR(b)*(x^2 - 2*xm*x), since the
C     ab-error covariance is COV(ab)=-xm*VAR(b)
C
C------------------------------------------------------------------------C
      DIMENSION X(NX),Y(NX),sig(NX)

      NITER = 2
      ITER = 0
      NTD = 0

10    sw=0.
      SX=0.
      SY=0.
      ST2=0.
      B=0.
      IF (iw.NE.0.0) THEN
        DO 11 I=1,NX
         if (sig(I).le.0. .or.x(i).lt.xmn.or.x(i).gt.xmx) goto 11
          wt = (1.0/sig(i)**2)
          sw=sw+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
        XM=SX/sw
        DO 12 I=1,NX
          if (sig(i).le.0. .or.x(i).lt.xmn.or.x(i).gt.xmx) goto 12
          Tx=(X(I)-XM)/sig(i)
          ST2=ST2+Tx*Tx
          B=B+Tx*Y(I)/sig(i)
12      CONTINUE
      else
        DO 13 I=1,NX
          if (x(i).lt.xmn .or. x(i).gt.xmx) goto 13
          sw = sw + 1.0
          SX=SX+X(I)
          SY=SY+Y(I)
13      CONTINUE
        XM=SX/sw
        DO 14 I=1,NX
          if (x(i).lt.xmn .or. x(i).gt.xmx) goto 14
          Tx=(X(I)-XM)
          ST2=ST2+Tx*Tx
          B=B+Tx*Y(I)
14      CONTINUE
      end if
      B=B/ST2
      A=(SY-SX*B)/sw
      SIGA=SQRT((1.+SX*SX/(sw*ST2))/sw)
      SIGB=SQRT(1./ST2)

C--   -- Estimate the Unbiased Standard-deviation.
      std = 0.0
      IF (iw.NE.0.0) THEN
        sw2 = 0.0
        DO I=1,NX
          if (sig(i).gt.0. .and. x(i).ge.xmn.and.x(i).le.xmx) then
            sw2 = sw2 + (1.0/sig(i)**4)
            std = std + ((y(i)-a-b*x(i))/sig(i))**2
          end if
        end do
        pts = sw**2/sw2
        std = sqrt(std/sw)
        if (pts.gt.2.5) std = sqrt(pts/(pts-2.0))*std
      ELSE
        do I=1,NX
         if (x(i).ge.xmn.and.x(i).le.xmx) std=std+(y(i)-a-b*x(i))**2
        end do
        pts = sw
        std = sqrt(std/(pts-2.0))
      END IF

C--   -- Eliminate highly deviant points and iterate
      if (iw.ne.0 .and.ITER.le.NITER .and.pts.gt.4. .and.T.gt.0) then
        ITER = ITER + 1
        nd = 0
        if (pts.le.3.0) then
         f = T*sqrt(sw/pts)
        else
         f = T*std
        end if
        do I=1,NX
         if (x(i).ge.xmn.and.x(i).le.xmx .and.sig(i).gt.0.
     &       .and. abs(y(i)-a-b*x(i)).gt.f) then
          sig(i) = -sig(i)
          nd = nd + 1
         end if
        end do
        ntd = ntd + nd
        if (nd.gt.0) goto 10
      end if

C--   -- Scale errors (assuming a good fit) if observed > expected error.
      if (pts.le.3.) then
       scale = 1.0
      else
       scale = std*sqrt(sw/pts)
      end if
      siga = siga*scale
      sigb = sigb*scale
      RETURN
      END

C------------------------------------------------------------------------C
      SUBROUTINE slp(X,Y,N,S,X0,Y0,DY0,XMN,XMX,A,B,DA,DB,CAB,STD,T,NTD)
C------------------------------------------------------------------------C
C
C     Fits the Slope b of line y = b*x + a, that passes through
C     X0,Y0 (whitin its error DY0), that minimizes 
C     b = sum{(y-y0)(x-x0)}/sum{(x-x0)^2},
C     We assume here that the given sigmas are a good estimator
C     of the real errors. We associate an error to the slope B
C     that is the larger of the one predicted by the observed
C     scatter (fit plus random error) and the one predicted
C     by the random errors (are the points alined by chance).
C
C     Dy0 - Sigma associated to Y0 (if any).
C
C     To exclude a Point, associate a non-possitive error to it.
C
C     THIS ROUTINE ITERATES ELIMINATING HIGHLY DEVIANT POINTS
C
C------------------------------------------------------------------------C
      implicit real*8 (A-H,O-Z)
      real*4 X(N),Y(N),S(N),X0,Y0,DY0,XMN,XMX,A,B,DA,DB,CAB,STD

      NITER = 2
      ITER = 0
      NTD = 0

10    SW = 0.
      SW2 = 0.
      SX = 0.
      SY = 0.
      SX2 = 0.
      SY2 = 0.
      SXY = 0.
      do i=1,N
        if (x(i).ge.xmn.and.x(i).le.xmx .and.S(i).gt.0.) then
          XI = (X(i)-X0)/S(i)
          YI = (Y(i)-Y0)/S(i)
          w = (1.0/S(i))**2
          SXY = SXY + YI*XI
          SX2 = SX2 + XI**2
          sw = sw + w
          sw2 = sw2 + w**2
          SX  = SX  + XI/S(i)
          SY  = SY  + YI/S(i)
          SY2 = SY2 + YI**2
        end if
      end do

C--   -- Optimaze the Ordinate at X0 within the errors of Y0.
c     dy = min(dy0,max(-dy0,(sx2*sy-sx*sxy)/(sw*sx2-sx**2)))
      dy = 0.0
      b  = (sxy-sx*dy)/sx2
      a  = Y0 + dy - b*X0
      db = sqrt(1.0/sx2)

!100   std = 0.0
      std = 0.0
      do I=1,N
        if (x(i).ge.xmn.and.x(i).le.xmx .and.S(i).gt.0.)
     &             std = std + ((y(i)-a-b*x(i))/S(i))**2
      end do
      pts = sw**2/sw2
      std = sqrt(std/sw)
      if (pts.gt.1.5) std = std*sqrt(pts/(pts-1.0))

C--   -- Eliminate deviators and iterate
      if (t.gt.0. .and. ITER.le.NITER .and. pts.gt.4.0) then
        ITER = ITER + 1
        nd = 0
        if (pts.le.3.0) then
         f = T*sqrt(sw/pts)
        else
         f = T*std
        end if
        do I=1,N
         if (x(i).ge.xmn.and.x(i).le.xmx .and.S(i).gt.0.
     &       .and. abs(y(i)-a-b*x(i)).gt.f) then
          S(i) = -S(i)
          nd = nd + 1
         end if
        end do
        ntd = ntd + nd
        if (nd.gt.0) goto 10
      end if

C--   -- Scale errors (assume a good fit) if observed error > expected.
C      scale = max(std*sqrt(sw/pts),1.0)
      if (pts.le.3.) then
       scale = 1.0
      else
       scale = std*sqrt(sw/pts)
      end if
      db = db*scale
      cab = -(x0*db**2)
      da = sqrt(dy0**2+(db*x0)**2)

      return
      end

