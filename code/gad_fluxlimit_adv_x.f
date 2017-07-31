C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_fluxlimit_adv_x.F,v 1.12 2008/02/29 030:59 mlosch Exp $

      SUBROUTINE GAD_FLUXLIMIT_ADV_X(iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskU,deltaT,deltaX,uVel,tracer,uT )
     
C !DESCRIPTION:
C Calculates the area integrated zonal flux due to advection of a tracer
C using second-order interpolation with a flux limiter:
C \begin{equation*}
C F^x_{adv} = U \overline{ \theta }^i
C - \frac{1}{2} \left(
C     [ 1 - \psi(C_r) ] |U|
C    + U \frac{u \Delta t}{\Delta x_c} \psi(C_r)
C              \right) \delta_i \theta
C \end{equation*}
C where the $\psi(C_r)$ is the limiter function and $C_r$ is
C the slope ratio.

      IMPLICIT NONE

      integer*4 i,j,k,t,im1,im2,ip1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaX(iMax,jMax)
      real*8 recipX(iMax,jMax)
      real*8 rA(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskU(iMax,jMax,Nr)
      real*8 uVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,uCFL,Cr,LCr
      real*8 uTrans,thetaMax,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(thetaMax = 1.D+20 )

      DO j=1,jMax
       DO i=1,iMax
          recipX(i,j)=1/deltaX(i,j)
       ENDDO
      ENDDO
      
      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax           
          IF (i .EQ. 2) THEN
           im1=1
           im2=iMax
          ELSEIF (i .EQ. 1) THEN
           im1=iMax
           im2=iMax-1
          ELSE
           im1=i-1
           im2=i-2
          ENDIF
          
          Rjp=(tracer(im1,j,k,t)-tracer( i ,j,k,t))
C     &     *maskU(im1,j,k)*maskC(im1,j,k)
          Rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
C     &     *maskU( i ,j,k)*maskC( i ,j,k)
          Rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
C     &     *maskU(im1,j,k)*maskC(im1,j,k)

          uCFL = ABS( uVel(i,j,k,t)*deltaT*recipX(i,j) )
          uTrans=uVel(i,j,k,t)*rA(i,j,k)

        IF (Rj.NE.0.) THEN
         IF (uTrans.GT.0) THEN
           Cr=Rjm/Rj
         ELSE
           Cr=Rjp/Rj
         ENDIF
        ELSE
         IF (uTrans.GT.0) THEN
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
        
C        uT(i,j,k,t) = maskU( i ,j,k)*maskU(im1,j,k)*
C     &                maskC( i ,j,k)*maskC(im1,j,k)* 
        uT(i,j,k,t) =   
     &   uTrans*(Tracer(i,j,k,t)+Tracer(im1,j,k,t))*half
     &   -ABS(uTrans)*((1.D0-LCr)+uCFL*LCr)
     &                    *Rj*half
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
