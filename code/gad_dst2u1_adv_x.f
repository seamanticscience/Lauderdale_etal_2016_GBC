C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst2u1_adv_x.F,v 1.8 2008/02/29 01:30:59 mlosch Exp $

      SUBROUTINE GAD_DST2U1_ADV_X(
     &           iMax,jMax,Nr,Nt,advectionScheme,rA,
     &           maskC,maskU,deltaT,deltaX,uVel,tracer,uT )

C !DESCRIPTION:
C  Calculates the area integrated zonal flux due to advection
C  of a tracer using second-order Direct Space and Time (DST-2)
C  interpolation (=Lax-Wendroff) or simple 1rst order upwind scheme.

      IMPLICIT NONE

      integer*4 advectionScheme
      integer*4 i,j,k,t,im1
      integer*4 iMax,jMax,Nr,Nt

      real*8 deltaT
      real*8 deltaX(iMax,jMax)
      real*8 recipX(iMax,jMax)
      real*8 rA    (iMax,jMax,Nr)
      real*8 maskC (iMax,jMax,Nr)
      real*8 maskU (iMax,jMax,Nr)
      real*8 uVel  (iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT    (iMax,jMax,Nr,Nt)

      real*8 uCFL, xLimit, uAbs, uTrans
CEOP

      IF ( advectionScheme.EQ.20 ) THEN
       xLimit = 1.D0
      ELSE
       xLimit = 0.D0
      ENDIF 

      DO j=1,jMax
       DO i=1,iMax
          recipX(i,j)=1/deltaX(i,j)
       ENDDO
      ENDDO
         
	  DO t=1,Nt		  
	   DO k=1,Nr		
	    DO j=1,jMax
         DO i=1,iMax
		  
		  IF (i .EQ. 1) THEN
		   im1=iMax
		  ELSE
		   im1=i-1
		  ENDIF
	  
		  uCFL = ABS( uVel(i,j,k,t)*deltaT*recipX(i,j) )

		  uTrans=uVel(i,j,k,t)*rA(i,j,k)
		
		  uAbs = ABS(uTrans)*( 1.D0 - xLimit*(1.D0 - uCFL) )
		  
		  uT(i,j,k,t) = 
     &		  ( uTrans+uAbs )* 0.5D0 * tracer(im1,j,k,t)
     &		+ ( uTrans-uAbs )* 0.5D0 * tracer(i  ,j,k,t)
	     ENDDO
	    ENDDO
	   ENDDO
	  ENDDO
      
      RETURN
      END
