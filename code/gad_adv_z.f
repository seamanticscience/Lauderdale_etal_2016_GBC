C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c2_adv_r.F,v 1.4 2001/09/21 13:143 adcroft Exp $

      SUBROUTINE GAD_ADV_Z( iMax,jMax,Nr,Nt,rA,maskC,
     &           wVel,tracer,wT )

C !DESCRIPTION:
C Calculates the area integrated vertical flux due to advection of a tracer
C using centered second-order interpolation:
C \begin{equation*}
C F^r_{adv} = W \overline{\theta}^k
C \end{equation*}

      IMPLICIT NONE

      integer*4 i,j,k,t,km1,kp1
      integer*4 iMax,jMax,Nr,Nt

      real*8 rA(iMax,jMax)
      real*8 maskC(iMax,jMax,Nr)
      real*8 wVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 wT(iMax,jMax,Nr,Nt)
C      real*8 af(iMax,jMax,Nr,Nt)
      real*8 half
      PARAMETER(half=1.D0/2.D0)
               
      DO t=1,Nt    
       DO k=1,Nr
        km1=max(1 ,k-1)
        kp1=min(Nr,k+1)
        
        IF ( k.EQ.1 .OR. k.GT.Nr) THEN
          DO j=1,jMax
           DO i=1,iMax
           wT(i,j,k,t) = 0.0000
          ENDDO
         ENDDO
        ELSE
          DO j=1,jMax
           DO i=1,iMax
           
C             wT(i,j,k,t) = maskC(i,j,km1)*maskC(i,j,k)*
             wT(i,j,k,t) = 
     &       wVel(i,j,k,t)*rA(i,j)
     &       *(Tracer(i,j,k,t)+Tracer(i,j,km1,t))*half
             
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      
      RETURN
      END
