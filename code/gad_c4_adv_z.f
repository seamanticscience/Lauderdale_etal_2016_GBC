C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_c4_adv_r.F,v 1.5 2003/05/01 16:12:36 jmc Exp $

      SUBROUTINE GAD_C4_ADV_Z( 
     I           iMax,jMax,Nr,Nt,
     I           rA, maskC,
     I           wVel, tracer,
     O           wT )

C !DESCRIPTION:
C Calculates the area integrated vertical flux due to advection of a tracer
C using centered fourth-order interpolation:
C \begin{equation*}
C F^r_{adv} = W \overline{ \theta - \frac{1}{6} \delta_{kk} \theta }^k
C \end{equation*}
C Near boundaries, the scheme reduces to a second if the flow is away
C from the boundary and to third order if the flow is towards
C the boundary.

      IMPLICIT NONE

      integer*4 i,j,k,t,km1,km2,kp1
      integer*4 iMax,jMax,Nr,Nt

      real*8 maskC(iMax,jMax,Nr)
      real*8 rA(iMax,jMax)
      real*8 wVel(iMax,jMax,Nr,Nt)
      real*8 tracer(iMax,jMax,Nr,Nt)
      real*8 wT(iMax,jMax,Nr,Nt)

      real*8 Rjm,Rj,Rjp,Rjjm,Rjjp
      real*8 rTrans,half,oneSixth,maskPM,maskBound
      PARAMETER(half=1.D0/2.D0)
      PARAMETER(oneSixth=1.D0/6.D0)

      DO t=1,Nt    
       DO k=1,Nr
        km2=MAX(1,k-2)
        km1=MAX(1,k-1)
        kp1=MIN(Nr,k+1)
        maskPM = 1.D0
        IF (k.LE.2 .OR. k.GE.Nr) maskPM = 0.D0
   
        IF ( k.EQ.1 .OR. k.GT.Nr) THEN
         DO j=1,jMax
          DO i=1,iMax
           wT(i,j,k,t) = 0.D0
          ENDDO
         ENDDO
        ELSE
         DO j=1,jMax
          DO i=1,iMax
           maskBound = maskPM*maskC(i,j,km2)*maskC(i,j,kp1)
           Rjp=(tracer(i,j,kp1,t)-tracer(i,j, k ,t))
C     &          *maskC(i,j,kp1)
           Rj =(tracer(i,j, k ,t)-tracer(i,j,km1,t))
           Rjm=(tracer(i,j,km1,t)-tracer(i,j,km2,t))
C     &          *maskC(i,j,km1)
           Rjjp=(Rjp-Rj)
           Rjjm=(Rj-Rjm)
           
           rTrans=wVel(i,j,k,t)*rA(i,j)

C           wT(i,j,k,t) = maskC(i,j,km1)*maskC(i,j,k)
C     &       *(
           wT(i,j,k,t) = 

     &       rTrans*(
     &         (Tracer(i,j,k,t)+Tracer(i,j,km1,t))*half
     &         -oneSixth*(Rjjm+Rjjp)*half 
     &              )
     &         +ABS(rTrans)*
     &         oneSixth*(Rjjm-Rjjp)*half*(1.D0 - maskBound)
C     &        )                      
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      
      RETURN
      END
