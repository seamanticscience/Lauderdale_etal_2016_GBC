C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_u3_adv_y.F,v 1.5 2011/03/29 15:50:30 jmc Exp $

      SUBROUTINE GAD_U3_ADV_Y( iMax,jMax,Nr,Nt,rA, 
     &              maskC,maskV,deltaT,deltaY,vVel,tracer,vT )

C !DESCRIPTION:
C Calculates the area integrated meridional flux due to advection of a tracer
C using upwind biased third-order interpolation (or the $\kappa=1/3$ scheme):
C \begin{equation*}
C F^y_{adv} = V \overline{ \theta  - \frac{1}{6} \delta_{jj} \theta }^j
C                 + \frac{1}{12} |V| \delta_{jjj} \theta
C \end{equation*}
C Near boundaries, mask all the gradients ==> still 3rd O.

      IMPLICIT NONE

      integer*4 i,j,k,t,jm1,jm2,jp1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaY(iMax,jMax)
      real*8 rA(iMax,jMax,Nr)
      real*8 maskV(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 vVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 vT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,Rjjm,Rjjp
      real*8 vTrans,vCFL,oneSixth,half
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(oneSixth=1.D0/6.D0)
      
      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          vT(i,j,k,t) = 0.d0
         ENDDO
        ENDDO
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
     &     *maskC(i,jp1,k)*maskC(i, j ,k)
          Rj =(tracer(i, j ,k,t)-tracer(i,jm1,k,t))
     &     *maskC(i, j ,k)*maskC(i,jm1,k)
          Rjm=(tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
     &     *maskC(i,jm1,k)*maskC(i,jm2,k)
          Rjjp=Rjp-Rj
          Rjjm=Rj-Rjm

          vTrans=vVel(i,j,k,t)*rA(i,j,k)
          
          vT(i,j,k,t) = maskC(i, j ,k)*
C     &                  maskC(i,jm1,k)*
C     &                  maskV(i, j ,k)*maskV(i,jm1,k)*     
C          vT(i,j,k,t) = 
     &   vTrans*(
     &     Tracer(i,j,k,t)+Tracer(i,jm1,k,t)-oneSixth*( Rjjp+Rjjm )
     &               )*half
     &  +ABS( vTrans )*half*oneSixth*( Rjjp-Rjjm )
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
