C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c2_adv_y.F,v 1.3 2001/09/21 13:143 adcroft Exp $

      SUBROUTINE GAD_ADV_Y( iMax,jMax,Nr,Nt,rA,maskC,maskV,
     &           vVel,tracer,vT )

C !DESCRIPTION:
C Calculates the area integrated meridional flux due to advection of a tracer
C using centered second-order interpolation:
C \begin{equation*}
C F^y_{adv} = V \overline{\theta}^j
C \end{equation*}

      IMPLICIT NONE

      integer*4 i,j,k,t,jm1,jp1
      integer*4 iMax,jMax,Nr,Nt
      real*8 maskC(iMax,jMax,Nr)
      real*8 maskV(iMax,jMax,Nr)
      real*8 rA(iMax,jMax,Nr)
      real*8 vVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)
      real*8 vT(iMax,jMax,Nr,Nt)
      real*8 half
      PARAMETER(half=1.D0/2.D0)      
      
         DO t=1,Nt
          DO k=1,Nr
           DO j=1,jMax
            DO i=1,iMax
             jm1=MAX(j-1, 1  )
             jp1=MIN(j+1,jMax)
             
             vT(i,j,k,t) = vVel(i,j,k,t)*rA(i,j,k)
     &             *(Tracer(i,j,k,t)+Tracer(i,jm1,k,t))*half

C             vT(i,j,k,t) = maskC(i,j,k)*maskC(i,jm1,k)*
C     &                     maskV(i,j,k)*maskV(i,jm1,k)*
C     &                     vVel(i,j,k,t)*rA(i,j,k)
C     &             *(Tracer(i,j,k,t)+Tracer(i,jm1,k,t))*half
            ENDDO
           ENDDO  
          ENDDO
         ENDDO

      RETURN
      END
