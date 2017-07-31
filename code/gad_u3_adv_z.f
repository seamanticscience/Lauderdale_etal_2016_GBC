C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_u3_adv_r.F,v 1.5 2014/08/18 12:22:46 jmc Exp $

      SUBROUTINE GAD_U3_ADV_Z( iMax,jMax,Nr,Nt,rA, 
     &              maskC,deltaT,deltaZ,wVel,tracer,wT )
     
C !DESCRIPTION:
C Calculates the area integrated vertical flux due to advection of a tracer
C using upwind biased third-order interpolation (or the $\kappa=1/3$ scheme):
C \begin{equation*}
C F^r_{adv} = W \overline{ \theta  - \frac{1}{6} \delta_{kk} \theta }^k
C                 + \frac{1}{12} |W| \delta_{kkk} \theta
C \end{equation*}
C Near boundaries, mask all the gradients ==> still 3rd O.

      IMPLICIT NONE

      integer*4 i,j,k,t,km1,km2,kp1
      integer*4 iMax,jMax,Nr,Nt
      
      real*8 deltaT
      real*8 deltaZ(Nr)
      real*8 rA(iMax,jMax)
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskUp(iMax,jMax,Nr)
      real*8 wVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 wT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,Rjjp,Rjjm
      real*8 rTrans,oneSixth,half
      PARAMETER(half=1.D0/2.D0)            
      PARAMETER(oneSixth=1.D0/6.D0)

      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax
          wT(i,j,k,t) = 0.d0
         ENDDO
        ENDDO
       ENDDO
      ENDDO             
      
      DO t=1,Nt    
       DO k=1,Nr
        km2=MAX(1,k-2)
        km1=MAX(1,k-1)
        kp1=MIN(Nr,k+1)

         DO j=1,jMax
          DO i=1,iMax
           Rjp=(tracer(i,j,kp1,t)-tracer(i,j, k ,t))
     &         *maskC(i,j,kp1)*maskC(i,j,k)
           Rj =(tracer(i,j, k ,t)-tracer(i,j,km1,t))
     &         *maskC(i,j,k)*maskC(i,j,km1)
           Rjm=(tracer(i,j,km1,t)-tracer(i,j,km2,t))
     &         *maskC(i,j,km1)*maskC(i,j,km2)
           
           Rjjp = Rjp-Rj
           Rjjm = Rj-Rjm

           rTrans= wVel(i,j,k,t)*rA(i,j)
         
           wT(i,j,k,t) = maskC(i,j,k)
C     &                  *maskC(i,j,km1)
C            wT(i,j,k,t) =
     &       *(rTrans*( (tracer(i,j,k,t)+tracer(i,j,km1,t))*half
     &                  -oneSixth*(Rjjm+Rjjp)*half  )
     &      +ABS(rTrans)*
     &                   oneSixth*(Rjjm-Rjjp)*half
     &                                  )
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
