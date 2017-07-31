C $Header: /u/gcmpack/MITgcm/pkg/generic_advdiff/gad_dst3_adv_r.F,v 1.11 2012/07/04 20:22:26 jmc Exp $

      SUBROUTINE GAD_DST3_ADV_Z( iMax,jMax,Nr,Nt,rA, 
     &              maskC,deltaT,deltaZ,wVel,tracer,wT )

C !DESCRIPTION:
C  Calculates the area integrated vertical flux due to advection of a tracer
C  using 3rd-order Direct Space and Time (DST-3) Advection Scheme

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

      real*8 Rjm,Rj,Rjp,d0,d1,oneSixth
      real*8 rTrans,wCFL,half
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
     &         *maskC(i,j,kp1)*maskC(i,j,k)
           Rj =(tracer(i,j,km1,t)-tracer(i,j,k,t))
     &         *maskC(i,j,k)*maskC(i,j,km1)
           Rjm=(tracer(i,j,km2,t)-tracer(i,j,km1,t))
     &         *maskC(i,j,km1)*maskC(i,j,km2)

           rTrans= wVel(i,j,k,t)*rA(i,j)
           wCFL=ABS( wVel(i,j,k,t)*deltaT*recipZ(k) )

           d0=(2.D0-wCFL)*(1.D0-wCFL)*oneSixth
           d1=(1.D0-wCFL*wCFL)*oneSixth

           wT(i,j,k,t)=maskC(i,j,k)*maskC(i,j,km1)           
C           wT(i,j,k,t)=(
     &    *(half*(rTrans+ABS(rTrans))
     &       *( tracer(i,j, k ,t) + (d0*Rj+d1*Rjp) )
     &   +half*(rTrans-ABS(rTrans))
     &       *( tracer(i,j,km1,t) - (d0*Rj+d1*Rjm) ))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      RETURN
      END
