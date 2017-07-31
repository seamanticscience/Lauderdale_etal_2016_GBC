C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_u3_adv_x.F,v 1.6 2011/03/29 15:50:30 jmc Exp $

      SUBROUTINE GAD_U3_ADV_X(iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskU,deltaT,deltaX,uVel,tracer,uT )

C !DESCRIPTION:
C Calculates the area integrated zonal flux due to advection of a tracer
C using upwind biased third-order interpolation (or the $\kappa=1/3$ scheme):
C \begin{equation*}
C F^x_{adv} = U \overline{ \theta  - \frac{1}{6} \delta_{ii} \theta }^i
C                 + \frac{1}{12} |U| \delta_{iii} \theta
C \end{equation*}
C Near boundaries, mask all the gradients ==> still 3rd O.

      IMPLICIT NONE

      integer*4 i,j,k,t,im1,im2,ip1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaX(iMax,jMax)
      real*8 rA(iMax,jMax,Nr)
      real*8 maskU(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 uVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 uT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,Rjjp,Rjjm
      real*8 uTrans,uCFL,oneSixth,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(oneSixth=1.D0/6.D0)

      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          uT(i,j,k,t) = 0.d0
         ENDDO
        ENDDO
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
          
          IF (i .EQ. iMax) THEN
           ip1=1
          ELSE
           ip1=i+1
          ENDIF          
          
          Rjp=(tracer(ip1,j,k,t)-tracer( i ,j,k,t))
     &     *maskC(ip1,j,k)*maskC( i ,j,k)
          Rj =(tracer( i ,j,k,t)-tracer(im1,j,k,t))
     &     *maskC( i ,j,k)*maskC(im1,j,k)
          Rjm=(tracer(im1,j,k,t)-tracer(im2,j,k,t))
     &     *maskC(im1,j,k)*maskC(im2,j,k)

          Rjjp=Rjp-Rj
          Rjjm=Rj-Rjm
        
          uTrans=uVel(i,j,k,t)*rA(i,j,k)
        
          uT(i,j,k,t) = maskC( i ,j,k)*
C     &                  maskC(im1,j,k)*
C     &                  maskU( i ,j,k)*maskU(im1,j,k)*     
C          uT(i,j,k,t) =
     &   uTrans*(
     &     Tracer(i,j,k,t)+Tracer(im1,j,k,t)-oneSixth*( Rjjp+Rjjm )
     &               )*half
     &     +ABS( uTrans )*half*oneSixth*( Rjjp-Rjjm )
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
