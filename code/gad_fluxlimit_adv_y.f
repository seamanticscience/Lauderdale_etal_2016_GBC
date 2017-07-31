C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_fluxlimit_adv_y.F,v 1.12 2007/04/04 039:06 jmc Exp $

      SUBROUTINE GAD_FLUXLIMIT_ADV_Y(iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskV,deltaT,deltaY,vVel,tracer,vT )

C !DESCRIPTION:
C Calculates the area integrated meridional flux due to advection of a tracer
C using second-order interpolation with a flux limiter:
C \begin{equation*}
C F^y_{adv} = V \overline{ \theta }^j
C - \frac{1}{2} \left(
C     [ 1 - \psi(C_r) ] |V|
C    + V \frac{v \Delta t}{\Delta y_c} \psi(C_r)
C              \right) \delta_j \theta
C \end{equation*}
C where the $\psi(C_r)$ is the limiter function and $C_r$ is
C the slope ratio.

      IMPLICIT NONE
      
      integer*4 i,j,k,t,jm1,jm2,jp1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaY(iMax,jMax)
      real*8 recipY(iMax,jMax)
      real*8 rA(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskV(iMax,jMax,Nr)
      real*8 vVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 vT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,vCFL,Cr,LCr
      real*8 vTrans,thetaMax,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(thetaMax = 1.D+20 )

      DO j=1,jMax
       DO i=1,iMax
          recipY(i,j)=1/deltaY(i,j)
       ENDDO
      ENDDO
      
      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax           
          jp1=MIN(j+1,jMax)
          jm1=MAX(j-1, 1  )
          jm2=MAX(j-2, 1  )
          
          Rjp=(tracer(i,jp1,k,t)-tracer(i, j ,k,t))
C     &     *maskV(i,jp1,k)*maskC(i,jp1,k)
          Rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
C     &     *maskV(i, j ,k)*maskC(i, j ,k)
          Rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
C     &     *maskV(i,jm1,k)*maskC(i,jm1,k)

          vCFL = ABS( vVel(i,j,k,t)*deltaT*recipY(i,j) )
          vTrans=vVel(i,j,k,t)*rA(i,j,k)

          IF (Rj.NE.0.) THEN
           IF (vTrans.GT.0) THEN
            Cr=Rjm/Rj
          ELSE
            Cr=Rjp/Rj
          ENDIF
         ELSE
          IF (vTrans.GT.0) THEN
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
             
C         vT(i,j,k,t) =maskV(i, j ,k)*maskV(i,jm1,k)*
C     &                maskC(i, j ,k)*maskC(i,jm1,k)* 
         vT(i,j,k,t) =   
     &   vTrans*(Tracer(i,j,k,t)+Tracer(i,jm1,k,t))*half
     &   -ABS(vTrans)*((1.D0-LCr)+vCFL*LCr)
     &                    *Rj*half
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
