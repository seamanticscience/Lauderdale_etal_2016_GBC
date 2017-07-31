C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst2u1_adv_r.F,v 1.4 2006/06/19 14:40:43 jmc Exp $

      SUBROUTINE GAD_DST2U1_ADV_Z(
     I           iMax,jMax,Nr,Nt,advectionScheme,rA,
     I           maskC,deltaT,deltaZ,wVel,tracer,wT )

C !DESCRIPTION:
C  Calculates the area integrated vertical flux due to advection
C  of a tracer using second-order Direct Space and Time (DST-2)
C  interpolation (=Lax-Wendroff) or simple 1rst order upwind scheme.

      IMPLICIT NONE

      integer*4 i,j,k,t,km1
      integer*4 iMax,jMax,Nr,Nt
      integer*4 advectionScheme

      real*8 deltaT
      real*8 deltaZ(Nr)
      real*8 recipZ(Nr)
      real*8 rA(iMax,jMax)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskUp(iMax,jMax,Nr)
      real*8 wVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 wT(iMax,jMax,Nr,Nt)
      
      real*8 wCFL, rLimit, wAbs, rTrans
CEOP
       real*8 rkSign
       PARAMETER(rkSign=-1.D0)
       
      IF ( advectionScheme.EQ.20 ) THEN
       rLimit = 1.D0
      ELSE
       rLimit = 0.D0
      ENDIF  

      DO k=1,Nr
       recipZ(k)=1/deltaZ(k)
      ENDDO
      
      DO t=1,Nt    
       DO k=1,Nr
        
        km1=MAX(1,k-1)
        
        IF ( k.LE.1 .OR. k.GT.Nr) THEN
         DO j=1,jMax
          DO i=1,iMax
           wT(i,j,k,t) = 0.D0
          ENDDO
         ENDDO
        ELSE
         DO j=1,jMax
          DO i=1,iMax
          wCFL = ABS( wVel(i,j,k,t)*deltaT*recipZ(k) )
          rTrans= wVel(i,j,k,t)*rA(i,j)

          wAbs = ABS(rTrans)*rkSign
     &       *( 1.D0 - rLimit*(1.D0 - wCFL) )
          wT(i,j,k,t) = maskC(i,j,km1)*(
     &             ( rTrans+wAbs )* 0.5D0 * tracer(i,j,km1,t)
     &           + ( rTrans-wAbs )* 0.5D0 * tracer(i,j,k  ,t)
     &                                  )
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END
