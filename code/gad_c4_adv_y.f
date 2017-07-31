C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c4_adv_y.F,v 1.6 2011/03/29 15:50:30 jmc Exp $

      SUBROUTINE GAD_C4_ADV_Y(
     I           iMax,jMax,Nr,Nt,
     I           rA, maskC, maskV,
     I           vVel,tracer,
     O           vT)

C !DESCRIPTION:
C Calculates the area integrated meridional flux due to advection of a tracer
C using centered fourth-order interpolation:
C \begin{equation*}
C F^y_{adv} = V \overline{ \theta - \frac{1}{6} \delta_{jj} \theta }^j
C \end{equation*}
C Near boundaries, the scheme reduces to a second if the flow is away
C from the boundary and to third order if the flow is towards
C the boundary.

      IMPLICIT NONE
      
      integer*4 i,j,k,t,jm1,jm2,jp1
      integer*4 iMax,jMax,Nr,Nt

      real*8 maskV(iMax,jMax,Nr)
      real*8 maskC(iMax,jMax,Nr)
      real*8 rA(iMax,jMax,Nr)
      real*8 vVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 vT(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,Rjjm,Rjjp
      real*8 vTrans,half,oneSixth
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(oneSixth=1.D0/6.D0)
      
      DO t=1,Nt    
       DO k=1,Nr
        DO j=1,jMax
         DO i=1,iMax           
          jp1=MIN(j+1,jMax)
          jm1=MAX(j-1, 1  )
          jm2=MAX(j-2, 1  )
       
          Rjp = (tracer(i,jp1,k,t)-tracer(i, j ,k,t))
C     &           *maskV(i,jp1,k)*maskC(i,jp1,k)
          Rj  = (tracer(i, j ,k,t)-tracer(i,jm1,k,t))
C     &           *maskV(i, j ,k)*maskC(i, j ,k)
          Rjm = (tracer(i,jm1,k,t)-tracer(i,jm2,k,t))
C     &           *maskV(i,jm1,k)*maskC(i,jm1,k)
          Rjjp=(Rjp-Rj)
          Rjjm=(Rj-Rjm)
          
          vTrans=vVel(i,j,k,t)*rA(i,j,k)

C          vT(i,j,k,t) = maskC(i, j ,k)*maskC(i,jm1,k)*
C     &                 *maskV(i, j ,k)*maskV(i,jm1,k)*     
          vT(i,j,k,t) = 
     &     vTrans*(
     &       Tracer(i,j,k,t)+Tracer(i,jm1,k,t)-oneSixth*( Rjjp+Rjjm )
     &            )*half
     &       +ABS( vTrans )*half*oneSixth*( Rjjp-Rjjm )
     &       *( 1.D0 - maskV(i,jm1,k)*maskV(i,jp1,k) )
         ENDDO
        ENDDO
       ENDDO
      ENDDO  

      RETURN
      END
