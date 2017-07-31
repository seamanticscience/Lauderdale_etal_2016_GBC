C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_fluxlimit_adv_r.F,v 1.12 2006/12/05 22:25:41 jmc Exp $

      SUBROUTINE GAD_FLUXLIMIT_ADV_Z( iMax,jMax,Nr,Nt,rA, 
     &              maskC,deltaT,deltaZ,wVel,tracer,wT )

C !DESCRIPTION:
C Calculates the area integrated vertical flux due to advection of a tracer
C using second-order interpolation with a flux limiter:
C \begin{equation*}
C F^x_{adv} = W \overline{ \theta }^k
C - \frac{1}{2} \left(
C     [ 1 - \psi(C_r) ] |W|
C    + W \frac{w \Delta t}{\Delta r_c} \psi(C_r)
C              \right) \delta_k \theta
C \end{equation*}
C where the $\psi(C_r)$ is the limiter function and $C_r$ is
C the slope ratio.

      IMPLICIT NONE

      integer*4 i,j,k,t,km1,km2,kp1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaZ(Nr)
      real*8 recipZ(Nr)
      real*8 rA(iMax,jMax)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskUp(iMax,jMax,Nr)
      real*8 wVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 wT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,wCFL,Cr,LCr
      real*8 rTrans,thetaMax,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(thetaMax = 1.D+20 )

      km2=MAX(1,k-2)
      km1=MAX(1,k-1)
      kp1=MIN(Nr,k+1)

      DO k=1,Nr
       recipZ(k)=1/deltaZ(k)
      ENDDO
      
      DO t=1,Nt    
       DO k=1,Nr
        km2=MAX(1,k-2)
        km1=MAX(1,k-1)
        kp1=MIN(Nr,k+1)

         DO j=1,jMax
          DO i=1,iMax
           Rjp=(tracer(i,j,k,t)-tracer(i,j,kp1,t))
C     &         *maskC(i,j,kp1)
           Rj =(tracer(i,j,km1,t)-tracer(i,j,k,t))
C     &         *maskC(i,j,k)*maskC(i,j,km1)
           Rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
C     &         *maskC(i,j,km1)

           rTrans= wVel(i,j,k,t)*rA(i,j)
		   wCFL = ABS(wVel(i,j,k,t)*deltaT*recipZ(k))

         IF (Rj.NE.0.) THEN
          IF (rTrans.LT.0.) THEN
            Cr=Rjm/Rj
          ELSE
            Cr=Rjp/Rj
          ENDIF
         ELSE
          IF (rTrans.LT.0.) THEN
            Cr=Rjm*thetaMax
          ELSE
            Cr=Rjp*thetaMax
          ENDIF
         ENDIF

C Statement function to describe flux limiter
C Upwind        Limiter(Cr)=0.
C Lax-Wendroff  Limiter(Cr)=1.
C Min-Mod       Limiter(Cr)=max(0.,min(1.,Cr))
C Suberbee      Limiter(Cr)=max(0.,max(min(1.,2*Cr),min(2.,Cr)))
	
c     Limiter(Cr)=0.
c     Limiter(Cr)=1.
c     Limiter(Cr)=max(0.D0,min(1.D0,Cr))
         LCr=max(0.D0,max(min(1.D0,2.D0*Cr),
     &                         min(2.D0,Cr)))
              
C         wT(i,j,k,t) = maskC(i,j,km1)*maskC(i,j,k)*(
         wT(i,j,k,t) = (
     &     rTrans*
     &        (tracer(i,j,k,t)+tracer(i,j,km1,t))*half
     &    +ABS(rTrans)*((1.D0-LCr)+wCFL*LCr)
     &                     *Rj*half )
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
