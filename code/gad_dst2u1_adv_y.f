C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst2u1_adv_y.F,v 1.7 2007/04/04 01:39:06 jmc Exp $

      SUBROUTINE GAD_DST2U1_ADV_Y(
     I           iMax,jMax,Nr,Nt,advectionScheme,rA,
     I           maskC,maskV,deltaT,deltaY,vVel,tracer,vT)

C !DESCRIPTION:
C  Calculates the area integrated meridional flux due to advection
C  of a tracer using second-order Direct Space and Time (DST-2)
C  interpolation (=Lax-Wendroff) or simple 1rst order upwind scheme.

      IMPLICIT NONE

      integer*4 advectionScheme
      integer*4 i,j,k,t,jm1
      integer*4 iMax,jMax,Nr,Nt

      real*8 deltaT
      real*8 deltaY(iMax,jMax)
      real*8 recipY(iMax,jMax)
      real*8 rA    (iMax,jMax,Nr)
      real*8 maskC (iMax,jMax,Nr)
      real*8 maskV (iMax,jMax,Nr)
      real*8 vVel  (iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 vT    (iMax,jMax,Nr,Nt)
      
      real*8 vCFL, yLimit, vAbs, vTrans
CEOP

      IF ( advectionScheme.EQ.20 ) THEN
       yLimit = 1.D0
      ELSE
       yLimit = 0.D0
      ENDIF
    
      DO j=1,jMax
       DO i=1,iMax
          recipY(i,j)=1/deltaY(i,j)
       ENDDO
      ENDDO
      
      DO t=1,Nt
       DO k=1,Nr
        DO j=1,jMax
         
         jm1=MAX(j-1, 1  )
         
         DO i=1,iMax
          
          vCFL = ABS( vVel(i,j,k,t)*deltaT*recipY(i,j) )
          vTrans=vVel(i,j,k,t)*rA(i,j,k)

          vAbs = ABS(vTrans)*( 1.D0 - yLimit*(1.D0 - vCFL) )
          vT(i,j,k,t) = ( vTrans+vAbs )* 0.5D0 * tracer(i,jm1,k,t)
     &                + ( vTrans-vAbs )* 0.5D0 * tracer(i,j  ,k,t)

         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      RETURN
      END
